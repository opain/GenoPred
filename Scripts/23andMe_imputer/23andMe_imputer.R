#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
# Large parts of were taken from ImputeMe (developed by Lasse Folkersen)
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--FID", action="store", default='FAM001', type='character',
		help="Family ID [optional]"),
make_option("--IID", action="store", default='ID001', type='character',
		help="Individual ID [optional]"),
make_option("--PID", action="store", default='0', type='character',
		help="Paternal ID [optional]"),
make_option("--MID", action="store", default='0', type='character',
		help="Maternal ID [optional]"),
make_option("--Pheno", action="store", default=-9, type='character',
		help="Phenotype: must be numeric. 1=control, 2=case [optional]"),
make_option("--Sex", action="store", default='i', type='character',
		help="Sex: i=infer from genotype, 1=male, 2=female, 0=missing [optional]"),
make_option("--geno", action="store", default=NULL, type='character',
		help="Path to 23andMe genetic data [required]"),
make_option("--plink", action="store", default='plink', type='character',
		help="Path to plink1.9 [optional]"),
make_option("--ref", action="store", default=NULL, type='character',
		help="Path to folder containing IMPUTE2 1KG reference data [required]"),
make_option("--shapeit", action="store", default='shapeit', type='character',
		help="Path to shapeit [optional]"),
make_option("--impute2", action="store", default='impute2', type='character',
		help="Path to impute2 [optional]"),
make_option("--output", action="store", default=NULL, type='character',
		help="Name of output files [required]"),
make_option("--n_core", action="store", default=1, type='numeric',
		help="Number of cores for parallel computing [optional]"),
make_option("--qctool", action="store", default='qctool', type='character',
		help="path to qctool2 [optional]"),
make_option("--chr", action="store", default=NULL, type='numeric',
		help="Chromosome number [optional]"),
make_option("--hardcall_thresh", action="store", default=0.9, type='numeric',
		help="path to qctool2 [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../Scripts/functions/misc.R')
library(foreach)
library(doMC)
registerDoMC(opt$n_core)

# Check required inputs
if(is.null(opt$geno)){
  stop('--geno must be specified.\n')
}
if(is.null(opt$ref)){
  stop('--ref must be specified.\n')
}
if(is.null(opt$output)){
  stop('--output must be specified.\n')
}
if(is.null(opt$qctool)){
  stop('--qctool must be specified.\n')
}
if(is.null(opt$chr)){
  stop('--chr must be specified.\n')
}

# Create output directory
opt$output_dir <- paste0(dirname(opt$output), '/')
system(paste0('mkdir -p ', opt$output_dir))

# Create temp directory
tmp_dir <- tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = '23andMe_imputer.R', start.time = start.time)

###
# Process input data
###

# If compressed, decompress
if(substr(opt$geno, (nchar(opt$geno) + 1) - 4, nchar(opt$geno)) == '.zip'){
	log <- system(paste0("unzip -d ", tmp_dir, " ", opt$geno), intern = T)
	opt$geno <- gsub(' ', '', gsub('.*inflating: ', '', log[grepl('inflating', log)]))
	log_add(log_file = log_file, message = 'Decompressed 23andMe genotype data.')
}

# Convert to PLINK format using plink
system(paste0(opt$plink,' --23file ',opt$geno,' ', opt$FID,' ',opt$IID,' ',opt$Sex,' ',opt$Pheno,' ',opt$PID,' ',opt$MID,' --out ',tmp_dir,'/geno --recode --chr ', opt$chr,' --memory 2000'))

# Read in map file
map <- fread(paste0(tmp_dir,'/geno.map'))
log_add(log_file = log_file, message = paste0('Input target data contains ', nrow(map), ' variants on chromosome ', opt$chr,'.'))

# Create list of SNPs that are duplicates
exclude <- map[duplicated(map[, 4]), 2]
log_add(log_file = log_file, message = paste0('Removing ',length(exclude),' SNPs that are duplicated.'))

write.table(exclude, file=paste0(tmp_dir,'/duplicates.snplist'), sep='\t', row.names=F, col.names=F, quote=F)
system(paste0(opt$plink,' --file ', tmp_dir, '/geno --recode --out ', tmp_dir, '/geno_nodup --exclude ', tmp_dir, '/duplicates.snplist'))

####
# Prepare for imputation using shapeit
####

log_add(log_file = log_file, message = 'Harmonsing and phasing genetic data using shapeit.')

# Check for strand flips
# We a introduce a dummy indvidual who is a heterozygote to make shapeit knows A2.
system(paste0(opt$shapeit, ' -check --input-ped ', tmp_dir, '/geno_nodup.ped ', tmp_dir, '/geno_nodup.map -M ', opt$ref, '/genetic_map_chr', opt$chr, '_combined_b37.txt --input-ref ',opt$ref,'/1000GP_Phase3_chr', opt$chr, '.hap.gz ', opt$ref,'/1000GP_Phase3_chr', opt$chr, '.legend.gz ', opt$ref, '/1000GP_Phase3.sample --output-log ', tmp_dir, '/geno_nodup.shapeit_log'))

# Read in the shapeit log file
logFile <- data.frame(fread(paste0(tmp_dir, '/geno_nodup.shapeit_log.snp.strand'), col.names = c("type1", "type2", "pos", "main_id", "main_A", "main_B", "main_flippable", "ref_id", "ref_A", "ref_B", "ref_flippable")))

# List SNPs that are missing, have mismatched alleles, or have missing alleles
omitMissing <- logFile$pos[logFile$type1 %in% 'Missing']
logStrand <- logFile[logFile$type1 %in% 'Strand',]
omitNonIdentical <- logStrand$pos[logStrand$main_A != logStrand$main_B]
omitBlank <- logStrand$pos[logStrand$main_A == '']

# Extract homozygote SNPs and remove INDELS
forceHomozygoteTable<-logStrand[
	logStrand$main_A == logStrand$main_B & 
	nchar(logStrand$ref_A) == 1 & 
	nchar(logStrand$ref_B) == 1 &
	!logStrand$main_A %in% c("D", "I") &
	!logStrand$main_B %in% c("D", "I"),]

# Remove SNPS where there are more than two alleles		
forceHomozygoteTable <- forceHomozygoteTable[sapply(apply(forceHomozygoteTable[, c('main_A', 'main_B', 'ref_A', 'ref_B')], 1, unique), length) == 2,]

# Remove variants with the same RSID
forceHomozygoteTable <- forceHomozygoteTable[!duplicated(forceHomozygoteTable[, 4]),]

# Read in the targrt sample map file
map <- data.frame(fread(paste0(tmp_dir, '/geno_nodup.map')))

# This adds a dummy individual who is heterozygous at every site
ped2 <- ped1 <- strsplit(readLines(paste0(tmp_dir, '/geno_nodup.ped'))," ")[[1]]
ped2[1] <- "Temporary"
ped2[2] <- "Non_person"
if((length(ped1) - 6) / 2 != dim(map)[1]){
	stop("mismatch between map and ped")
}
replacementPos <- which(map$V2 %in% forceHomozygoteTable$main_id)
A1_pos <- 7 + 2 * (replacementPos - 1)
A2_pos <- 8 + 2 * (replacementPos - 1)
ped2[A1_pos] <- forceHomozygoteTable[, 9]
ped2[A2_pos] <- forceHomozygoteTable[, 10]
ped <- rbind(ped1, ped2)
write.table(ped, paste0(tmp_dir, '/geno_nodup.ped'), sep=" ", col.names=F, row.names=F, quote=F)
omitRemaining <- logStrand[!logStrand[, 4] %in% forceHomozygoteTable[, 4], 3]
write.table(c(omitNonIdentical, omitBlank, omitMissing, omitRemaining), file = paste0(tmp_dir, '/exclusions.txt'), sep = '\t', row.names = F, col.names = F, quote = F)
	
# Run shapeit command
system(paste0(opt$shapeit, ' --input-ped ', tmp_dir, '/geno_nodup.ped ', tmp_dir, '/geno_nodup.map -M ', opt$ref,'/genetic_map_chr', opt$chr ,'_combined_b37.txt --input-ref ', opt$ref, '/1000GP_Phase3_chr', opt$chr, '.hap.gz ', opt$ref, '/1000GP_Phase3_chr', opt$chr, '.legend.gz ', opt$ref, '/1000GP_Phase3.sample --output-log ',tmp_dir, '/geno_nodup.shapeit_log --exclude-snp ',tmp_dir, '/exclusions.txt -O ',tmp_dir, '/geno_nodup.harmonised_temp'))

# Remove dummy individual
system(paste0("cut --delimiter=' ' -f 1-7 ",tmp_dir, '/geno_nodup.harmonised_temp.haps > ',tmp_dir, '/geno_nodup.harmonised.haps'))
system(paste0("head -n 3 ",tmp_dir, '/geno_nodup.harmonised_temp.sample > ',tmp_dir, '/geno_nodup.harmonised.sample'))

###########################################
# Impute genetic data
###########################################

log_add(log_file = log_file, message = 'Imputing data using Impute2.')

# Look up length of chromosome
maxPos <- system(paste0('zcat ', opt$ref, '/1000GP_Phase3_chr', opt$chr, ".legend.gz | tail -n 1 | cut -d ' ' -f 2"),intern=T)

# Create a data.frame containing list all 5Mb chunks for imputation
starts <- seq(0, as.numeric(maxPos), 5e6)
ends <- starts + 5e6

foreach(i=1:length(starts), .combine=c, .options.multicore=list(preschedule=FALSE)) %dopar% {
	imp_log <- system(paste0(opt$impute2,' -m ', opt$ref, '/genetic_map_chr', opt$chr, '_combined_b37.txt -h ', opt$ref, '/1000GP_Phase3_chr', opt$chr, '.hap.gz -l ', opt$ref, '/1000GP_Phase3_chr', opt$chr, '.legend.gz -known_haps_g ', tmp_dir, '/geno_nodup.harmonised.haps -int ', starts[i], ' ', ends[i], ' -Ne 20000 -o ', tmp_dir, '/geno_nodup.', starts[i], '_', ends[i]))
    
	# test for memory-lack bug (error will be 137 if killed, otherwise 0)
    if(imp_log == 137){
		log_add(log_file = log_file, message = paste0("Chunk ", starts[i],'-',ends[i],": There was a memory problem."))
		stop('There was a memory problem.')
	}
}

# Combine the chunks
chunk_files <- paste0(tmp_dir, '/geno_nodup.', starts, '_', ends)

# Check which chunk_files exist
chunk_files_present<-NULL
for(k in 1:length(chunk_files)){
	if(file.exists(chunk_files[k])){
		chunk_files_present<-c(chunk_files_present,chunk_files[k])
	}
}

system(paste0('cat ', paste(chunk_files_present, collapse=' '),' > ',tmp_dir, '/geno_nodup_imp.gen'))

# Move log files to a log file folder
system(paste0('mkdir -p ', opt$output_dir, '/impute2_logs/chr', opt$chr))
for(k in 1:length(chunk_files)){
	system(paste0('mv ',chunk_files[k],'_summary ',opt$output_dir, '/impute2_logs/chr', opt$chr,'/'))
	system(paste0('mv ',chunk_files[k],'_warnings ',opt$output_dir, '/impute2_logs/chr', opt$chr,'/'))
}

#####
# Convert into PLINK format
#####

log_add(log_file = log_file, message = 'Converting to PLINK format.')

# Insert chromosome number
gen<-fread(paste0(tmp_dir, '/geno_nodup_imp.gen'))
gen$V1<-as.character(opt$chr)

# Exclude INDELS
gen<-gen[nchar(gen$V4) == 1 & nchar(gen$V5) == 1,]

# Include only the RSID if available
gen$V2[grepl('rs',gen$V2)]<-gsub(':.*','',gen$V2[grepl('rs',gen$V2)])
fwrite(gen, paste0(tmp_dir, '/geno_nodup_imp.gen'), col.names=F, row.names=F, quote=F, sep=' ')

# Save as a plink file
# Make sure qctool will find the library
system('export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}')
system(paste0(opt$qctool,' -g ',tmp_dir, '/geno_nodup_imp.gen -threshold ',opt$hardcall_thresh,' -og ', tmp_dir, '/geno_nodup_imp -ofiletype binary_ped'))

# Remove missing SNPs in the plink files
system(paste0(opt$plink, ' --bfile ', tmp_dir, '/geno_nodup_imp --geno 0.1 --make-bed --out ', opt$output, ' --memory 2000'))

# Insert the chromosome number to the bim file
bim<-fread(paste0(opt$output,'.bim'))
bim$V1 <- opt$chr
fwrite(bim, paste0(opt$output,'.bim'), col.names=F, row.names=F, quote=F, sep='\t')
	
# Rename the gen file
system(paste0('mv ',tmp_dir, '/geno_nodup_imp.gen ',opt$output,'.gen'))

# Update fam file with IDs and inferred sex
system(paste0("cut -f 1-6 -d ' ' ",tmp_dir, '/geno.ped > ', opt$output,'.fam'))

# Create a sample file for the bgen files
sample_file<-data.frame(ID=c(0,paste0(opt$FID,'_',opt$IID)),
						sex=c('D','M'))			

write.table(sample_file, paste0(opt$output,'.sample'), col.names=T, row.names=F, quote=F)

# Count the number of variants after imputation in gen files and plink files
final_n_gen<-gsub(' .*', '', system(paste0('wc -l ',opt$output,'.chr*.gen'), intern=T))
final_n_plink<-gsub(' .*','', system(paste0('wc -l ',opt$output,'.chr*.bim'), intern=T))

log_add(log_file = log_file, message = paste0('.gen files contain ',final_n_gen,' SNPs.'))
log_add(log_file = log_file, message = paste0('.plink files contain ',final_n_plink,' SNPs.'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()
