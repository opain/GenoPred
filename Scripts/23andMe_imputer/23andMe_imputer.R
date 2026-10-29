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
make_option("--plink2", action="store", default='plink2', type='character',
    help="Path to plink2 [optional]"),
make_option("--ref", action="store", default=NULL, type='character',
		help="Path to folder containing IMPUTE2-format 1KG reference data, used for phasing only [required]"),
make_option("--shapeit", action="store", default='shapeit', type='character',
		help="Path to shapeit [optional]"),
make_option("--bcftools", action="store", default='bcftools', type='character',
		help="Path to bcftools [optional]"),
make_option("--impute5", action="store", default='impute5', type='character',
		help="Path to impute5 [optional]"),
make_option("--ref_name", action="store", default='1000GP_Phase3', type='character',
		help="File stem of the reference panel in --ref and --ref5 (<stem>_chr<N>.hap.gz/.legend.gz, <stem>.sample, <stem>_chr<N>_xcf.bcf) [optional]"),
make_option("--ref5", action="store", default=NULL, type='character',
		help="Path to folder containing the IMPUTE5 xcf-format 1KG reference data (1000GP_Phase3_chr<N>_xcf.bcf), used for imputation [required]"),
make_option("--map5", action="store", default=NULL, type='character',
		help="Path to folder containing IMPUTE5/SHAPEIT4-format genetic maps (chr<N>.b37.gmap.gz), used for imputation [required]"),
make_option("--output", action="store", default=NULL, type='character',
		help="Name of output files [required]"),
make_option("--n_core", action="store", default=1, type='numeric',
		help="Number of cores for parallel computing [optional]"),
make_option("--chr", action="store", default=NULL, type='numeric',
		help="Chromosome number [optional]"),
make_option("--hardcall_thresh", action="store", default=0.9, type='numeric',
		help="Hardcall threshold [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')

# Check required inputs
if(is.null(opt$geno)){
  stop('--geno must be specified.\n')
}
if(is.null(opt$ref)){
  stop('--ref must be specified.\n')
}
if(is.null(opt$ref5)){
  stop('--ref5 must be specified.\n')
}
if(is.null(opt$map5)){
  stop('--map5 must be specified.\n')
}
if(is.null(opt$output)){
  stop('--output must be specified.\n')
}
if(is.null(opt$chr)){
  stop('--chr must be specified.\n')
}

# Create output directory
opt$output_dir <- paste0(dirname(opt$output), '/')
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# Create temp directory
tmp_dir <- tempdir()

run_tool <- function(command, args, stdout = '') {
  status <- system2(command, shQuote(as.character(args)), stdout = stdout)
  if(status != 0) stop(command, ' failed with exit status ', status)
}
tmp <- function(name) file.path(tmp_dir, name)

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = '23andMe_imputer.R', start.time = start.time)

###
# Process input data
###

# If compressed, decompress
if(substr(opt$geno, (nchar(opt$geno) + 1) - 4, nchar(opt$geno)) == '.zip'){
	entries <- utils::unzip(opt$geno, list = TRUE)$Name
	entries <- entries[!grepl('/$', entries)]
	if(length(entries) != 1 || grepl('^/|(^|/)\\.\\.(/|$)', entries[1])) stop('Genotype zip must contain one safe file path.')
	utils::unzip(opt$geno, files = entries, exdir = tmp_dir)
	opt$geno <- file.path(tmp_dir, entries[1])
	log_add(log_file = log_file, message = 'Decompressed 23andMe genotype data.')
}

# Convert to PLINK format using plink
run_tool(opt$plink, c('--23file', opt$geno, opt$FID, opt$IID, opt$Sex, opt$Pheno, opt$PID, opt$MID, '--out', tmp('geno'), '--recode', '--chr', opt$chr, '--memory', 2000))

# Read in map file
map <- fread(paste0(tmp_dir,'/geno.map'))
log_add(log_file = log_file, message = paste0('Input target data contains ', nrow(map), ' variants on chromosome ', opt$chr,'.'))

# Create list of SNPs that are duplicates
exclude <- map[duplicated(map[, 4]), 2]
log_add(log_file = log_file, message = paste0('Removing ',length(exclude),' SNPs that are duplicated.'))

write.table(exclude, file=paste0(tmp_dir,'/duplicates.snplist'), sep='\t', row.names=F, col.names=F, quote=F)
run_tool(opt$plink, c('--file', tmp('geno'), '--recode', '--out', tmp('geno_nodup'), '--exclude', tmp('duplicates.snplist')))

####
# Prepare for imputation using shapeit
####

log_add(log_file = log_file, message = 'Harmonsing and phasing genetic data using shapeit.')

# Check for strand flips
# We a introduce a dummy indvidual who is a heterozygote to make shapeit knows A2.
shapeit_input <- c('--input-ped', tmp('geno_nodup.ped'), tmp('geno_nodup.map'), '-M', file.path(opt$ref, paste0('genetic_map_chr', opt$chr, '_combined_b37.txt')), '--input-ref', file.path(opt$ref, paste0(opt$ref_name, '_chr', opt$chr, '.hap.gz')), file.path(opt$ref, paste0(opt$ref_name, '_chr', opt$chr, '.legend.gz')), file.path(opt$ref, paste0(opt$ref_name, '.sample')), '--output-log', tmp('geno_nodup.shapeit_log'))
run_tool(opt$shapeit, c('-check', shapeit_input))

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
run_tool(opt$shapeit, c(shapeit_input, '--exclude-snp', tmp('exclusions.txt'), '-O', tmp('geno_nodup.harmonised_temp')))

# Remove dummy individual
run_tool('cut', c('-d', ' ', '-f', '1-7', tmp('geno_nodup.harmonised_temp.haps')), tmp('geno_nodup.harmonised.haps'))
writeLines(readLines(tmp('geno_nodup.harmonised_temp.sample'), n=3), tmp('geno_nodup.harmonised.sample'))

###########################################
# Impute genetic data using IMPUTE5
###########################################

log_add(log_file = log_file, message = 'Imputing data using IMPUTE5.')

# IMPUTE5 needs the phased target as indexed VCF/BCF, not SHAPEIT's native
# .haps/.sample. bcftools' --hapsample2vcf requires the ID in column 1 or 2 to
# be "CHR:POS_REF_ALT" (underscore before the alleles - see `man bcftools`);
# SHAPEIT's own IDs are plain rsIDs, so column 2 is rebuilt from the chr/pos/
# allele columns already present in the .haps file.
run_tool('awk', c('{$2=$1":"$3"_"$4"_"$5; print}', tmp('geno_nodup.harmonised.haps')), tmp('geno_nodup.harmonised.haps_id'))
run_tool(opt$bcftools, c('convert', '--hapsample2vcf', paste0(tmp('geno_nodup.harmonised.haps_id'), ',', tmp('geno_nodup.harmonised.sample')), '-O', 'b', '-o', tmp('target.bcf')))

# xcftools/impute5 both require an INFO/AC field, which bcftools convert's
# output doesn't carry - same fix as the one-time reference conversion
# (hpc/setup/04_build_impute5_reference.sh).
run_tool(opt$bcftools, c('+fill-tags', tmp('target.bcf'), '-O', 'b', '-o', tmp('target_tagged.bcf'), '--', '-t', 'AC,AN'))
run_tool(opt$bcftools, c('index', '-f', tmp('target_tagged.bcf')))

# Whole-chromosome region for --r/--buffer-region (IMPUTE5 imputes the full
# chromosome in one call - no manual chunking needed, unlike IMPUTE2). Both
# --r and --buffer-region are required by impute5 even when there's no actual
# buffer beyond the chromosome itself, and need an explicit position range,
# not just the bare chromosome ID.
legend <- gzfile(file.path(opt$ref, paste0(opt$ref_name, '_chr', opt$chr, '.legend.gz')), 'rt')
last_line <- ''
repeat { lines <- readLines(legend, n = 100000); if(!length(lines)) break; last_line <- lines[length(lines)] }
close(legend)
maxPos <- strsplit(trimws(last_line), '[[:space:]]+')[[1]][2]
if(is.na(suppressWarnings(as.numeric(maxPos)))) stop('Invalid reference legend position.')
region <- paste0(opt$chr, ':1-', maxPos)

ref5_bcf <- paste0(opt$ref5, '/', opt$ref_name, '_chr', opt$chr, '_xcf.bcf')
map5_file <- paste0(opt$map5, '/chr', opt$chr, '.b37.gmap.gz')

run_tool(opt$impute5, c('--h', ref5_bcf, '--g', tmp('target_tagged.bcf'), '--o', tmp('imputed.bcf'), '--r', region, '--buffer-region', region, '--m', map5_file, '--threads', opt$n_core))

# Convert the imputed BCF back to the classic .gen format the rest of this
# script (and the pipeline downstream) expects. --tag GP extracts genotype
# probabilities, matching IMPUTE2's own .gen semantics.
run_tool(opt$bcftools, c('convert', '--gensample', tmp('geno_nodup_imp'), '--tag', 'GP', tmp('imputed.bcf')))
run_tool('gunzip', c('-f', tmp('geno_nodup_imp.gen.gz')))

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

# Create a sample file
write.table(
  data.frame(
    ID_1 = c(0, 'X'),
    ID_2 = c(0, 'X'),
    missing = c(0, NA),
    sex = c('D', NA)
  ),
  paste0(tmp_dir, '/geno_nodup_imp.sample'),
  col.names = T,
  row.names = F,
  quote = F
)

# Save as a plink file
run_tool(opt$plink2, c('--gen', tmp('geno_nodup_imp.gen'), 'ref-last', '--sample', tmp('geno_nodup_imp.sample'), '--hard-call-threshold', 0.1, '--make-bed', '--memory', 5000, '--threads', 1, '--out', tmp('geno_nodup_imp')))

# Remove missing SNPs in the plink files
run_tool(opt$plink, c('--bfile', tmp('geno_nodup_imp'), '--geno', 0.1, '--make-bed', '--out', tmp('geno_nodup_imp_geno'), '--memory', 2000))

# Insert the chromosome number to the bim file
bim<-fread(paste0(tmp_dir, '/geno_nodup_imp_geno.bim'))
bim$V1 <- opt$chr
fwrite(bim, paste0(tmp_dir, '/geno_nodup_imp_geno.bim'), col.names=F, row.names=F, quote=F, sep='\t')

# Rename the imputed genotype data to output
for(ext in c('bed', 'bim', 'gen')) {
  source <- tmp(paste0(if(ext == 'gen') 'geno_nodup_imp' else 'geno_nodup_imp_geno', '.', ext))
  destination <- paste0(opt$output, '.', ext)
  unlink(destination)
  if(!file.copy(source, destination)) stop('Could not copy output: ', destination)
  unlink(source)
}

# Update fam file with IDs and inferred sex
unlink(paste0(opt$output, '.fam'))
run_tool('cut', c('-f', '1-6', '-d', ' ', tmp('geno.ped')), paste0(opt$output, '.fam'))

# Create a sample file for the bgen files
sample_file<-data.frame(
  ID=c(0,paste0(opt$FID,'_',opt$IID)),
	sex=c('D', NA))

unlink(paste0(opt$output, '.sample'))
write.table(sample_file, paste0(opt$output,'.sample'), col.names=T, row.names=F, quote=F)

# Count the number of variants after imputation in gen files and plink files
final_n_gen <- strsplit(trimws(system2('wc', c('-l', shQuote(paste0(opt$output, '.gen'))), stdout=TRUE)), '[[:space:]]+')[[1]][1]
final_n_plink <- nrow(bim)

log_add(log_file = log_file, message = paste0('.gen files contain ',final_n_gen,' SNPs.'))
log_add(log_file = log_file, message = paste0('.plink files contain ',final_n_plink,' SNPs.'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()
