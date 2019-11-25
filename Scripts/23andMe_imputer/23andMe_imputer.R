#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

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
		help="Sex: 1=male, 2=female, 0=missing [optional]"),
make_option("--geno", action="store", default=NA, type='character',
		help="Path to 23andMe genetic data [required]"),
make_option("--plink", action="store", default=NA, type='character',
		help="Path to plink1.9 [required]"),
make_option("--ref", action="store", default=NA, type='character',
		help="Path to folder containing IMPUTE2 1KG reference data [required]"),
make_option("--shapeit", action="store", default=NA, type='character',
		help="Path to shapeit [required]"),
make_option("--impute2", action="store", default=NA, type='character',
		help="Path to impute2 [required]"),
make_option("--out", action="store", default=NA, type='character',
		help="Name of output files [required]"),
make_option("--n_core", action="store", default=1, type='numeric',
		help="Number of cores for parallel computing [required]"),
make_option("--qctool", action="store", default='NA', type='character',
		help="path to qctool2 [required]"),
make_option("--hardcall_thresh", action="store", default=0.9, type='numeric',
		help="path to qctool2 [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Create chromosome variable
CHROMS<-c(sample(c('X',1:22)))

sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = F)
cat(
"#################################################################
# 23andMe_imputer.R
# This is an adaptation of Lasse Folkerson's impute.me code
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at",as.character(start.time),'
Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

library(data.table)
library(foreach)
library(doMC)
registerDoMC(opt$n_core)

# If compressed, decompress
if(substr(opt$geno,(nchar(opt$geno)+1)-4,nchar(opt$geno)) == '.zip'){
  system(paste0("unzip -d ",dirname(opt$geno)," ",opt$geno))
  opt$geno<-sub('.zip','.txt',opt$geno)
  sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
	cat('Decompressed 23andMe genotype data.\n')
  sink()
}

# Convert to PLINK format using plink
system(paste0(opt$plink,' --23file ',opt$geno,' ', opt$FID,' ',opt$IID,' ',opt$Sex,' ',opt$Pheno,' ',opt$PID,' ',opt$MID,' --out ',opt$out,' --recode --memory 2000'))

# Count the number of SNPs
n_geno<-system(paste0('wc -l ',opt$out,'.map'),intern=T)
n_geno<-unlist(strsplit(n_geno,' '))[1]
sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
cat('Input target data contains',n_geno,'variants\n')
sink()

# Create list of SNPs that are duplicates
map<-data.frame(fread(paste0(opt$out,'.map')))
exclude<-map[duplicated(map[,4]),2]
sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
cat(paste0('Removed ',length(exclude),' SNPs that were duplicated.\n'))
sink()
write.table(exclude,file=paste0(opt$out,'_duplicates.snplist'),sep='\t',row.names=F,col.names=F,quote=F)

sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
cat('Harmonsing and phasing genetic data using shapeit...')
sink()

# Loop over each chromosomes
foreach(i=CHROMS, .combine=c) %dopar% {
	
	#First in loop - extract only one specific chromosome
	log1<-system(paste0(opt$plink,' --file ',opt$out,' --chr ',i,' --recode --out ',opt$out,'.chr',i,' --exclude ',opt$out,'_duplicates.snplist'))

	#if X chromosome is missing it is allowed to skip forward
	if(log1 == 13 & i == "X"){
		sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
		cat(paste0("Didn't find X-chr data, so skipping that.\n"))
		sink()
	}
	
	# Check for strand flips
	# We a introduce a dummy indvidual who is a heterozygote to make shapeit knows A2.
	if(i == 'X'){
		system(paste0(opt$shapeit,' -check --input-ped ',opt$out,'.chr',i,'.ped ',opt$out,'.chr',i,'.map -M ',opt$ref,'/genetic_map_chr',i,'_nonPAR_combined_b37.txt --input-ref ',opt$ref,'/1000GP_Phase3_chr',i,'_NONPAR.hap.gz ',opt$ref,'/1000GP_Phase3_chr',i,'_NONPAR.legend.gz ',opt$ref,'/1000GP_Phase3.sample --output-log ',opt$out,'.chr',i,'.shapeit_log'))
	} else {
		system(paste0(opt$shapeit,' -check --input-ped ',opt$out,'.chr',i,'.ped ',opt$out,'.chr',i,'.map -M ',opt$ref,'/genetic_map_chr',i,'_combined_b37.txt --input-ref ',opt$ref,'/1000GP_Phase3_chr',i,'.hap.gz ',opt$ref,'/1000GP_Phase3_chr',i,'.legend.gz ',opt$ref,'/1000GP_Phase3.sample --output-log ',opt$out,'.chr',i,'.shapeit_log'))
	}	
	# Read in the shapeit log file
	logFile<-data.frame(fread(paste0(opt$out,'.chr',i,'.shapeit_log.snp.strand'), col.names=c("type1","type2","pos","main_id","main_A","main_B","main_flippable","ref_id","ref_A","ref_B","ref_flippable")))
	
	# List SNPs that are missing, have mismatched alleles, or have missing alleles
	omitMissing<-logFile$pos[logFile$type1 %in% 'Missing']
	logStrand<-logFile[logFile$type1 %in% 'Strand',]
	omitNonIdentical<-logStrand$pos[logStrand$main_A != logStrand$main_B]
	omitBlank<-logStrand$pos[logStrand$main_A == '']
	
	# Extract homozygote SNPs and remove INDELS
	forceHomozygoteTable<-logStrand[
		logStrand$main_A == logStrand$main_B & 
			nchar(logStrand$ref_A)==1 & 
			nchar(logStrand$ref_B)==1 &
			!logStrand$main_A %in% c("D","I") &
			!logStrand$main_B %in% c("D","I") 
		,]
	
	# Remove SNPS where there are more than two alleles		
	forceHomozygoteTable<-forceHomozygoteTable[sapply(apply(forceHomozygoteTable[,c('main_A','main_B','ref_A','ref_B')],1,unique),length)==2,]
	
	# Remove variants with the same RSID
	forceHomozygoteTable<-forceHomozygoteTable[!duplicated(forceHomozygoteTable[,4]),]
	
	# Read in the targrt sample map file
	map<-data.frame(fread(paste0(opt$out,'.chr',i,'.map')))
	
	# This adds a dummy individual who is heterozygous at every site
	ped2<-ped1<-strsplit(readLines(paste0(opt$out,'.chr',i,'.ped'))," ")[[1]]
	ped2[1]<-"Temporary"
	ped2[2]<-"Non_person"
	if((length(ped1)-6) / 2 != dim(map)[1]){
			stop("mismatch between map and ped")
	}
	replacementPos<-which(map$V2 %in% forceHomozygoteTable$main_id)
	A1_pos<-7+2*(replacementPos-1)
	A2_pos<-8+2*(replacementPos-1)
	ped2[A1_pos]<-forceHomozygoteTable[,9]
	ped2[A2_pos]<-forceHomozygoteTable[,10]
	ped<-rbind(ped1,ped2)
	write.table(ped,paste0(opt$out,'.chr',i,'.ped',sep=""),sep=" ",col.names=F,row.names=F,quote=F)
	omitRemaining<-logStrand[!logStrand[,4]%in%forceHomozygoteTable[,4],3]
	write.table(c(omitNonIdentical,omitBlank,omitMissing,omitRemaining),file=paste0(opt$out,'.chr',i,'_exclusions'),sep='\t',row.names=F,col.names=F,quote=F)
	
	# Run shapeit command
	if(i == 'X'){
		system(paste0(opt$shapeit,' --input-ped ',opt$out,'.chr',i,'.ped ',opt$out,'.chr',i,'.map -M ',opt$ref,'/genetic_map_chr',i,'_nonPAR_combined_b37.txt --input-ref ',opt$ref,'/1000GP_Phase3_chr',i,'_NONPAR.hap.gz ',opt$ref,'/1000GP_Phase3_chr',i,'_NONPAR.legend.gz ',opt$ref,'/1000GP_Phase3.sample --output-log ',opt$out,'.chr',i,'.shapeit_log --exclude-snp ',opt$out,'.chr',i,'_exclusions -O ',opt$out,'.chr',i,'.harmonised_temp'))
	} else {
		system(paste0(opt$shapeit,' --input-ped ',opt$out,'.chr',i,'.ped ',opt$out,'.chr',i,'.map -M ',opt$ref,'/genetic_map_chr',i,'_combined_b37.txt --input-ref ',opt$ref,'/1000GP_Phase3_chr',i,'.hap.gz ',opt$ref,'/1000GP_Phase3_chr',i,'.legend.gz ',opt$ref,'/1000GP_Phase3.sample --output-log ',opt$out,'.chr',i,'.shapeit_log --exclude-snp ',opt$out,'.chr',i,'_exclusions -O ',opt$out,'.chr',i,'.harmonised_temp'))
	}
	
	# Remove dumy individual
	system(paste0("cut --delimiter=' ' -f 1-7 ",opt$out,'.chr',i,'.harmonised_temp.haps > ',opt$out,'.chr',i,'.harmonised.haps'))
	system(paste0("head -n 3 ",opt$out,'.chr',i,'.harmonised_temp.sample > ',opt$out,'.chr',i,'.harmonised.sample',sep=""))

	# Delete temporary files 
	system(paste0('rm ',opt$out, '.chr',i,'.shapeit_log.snp.strand.exclude'))
	system(paste0('rm ',opt$out, '.chr',i,'.shapeit_log.snp.strand'))
	system(paste0('rm ',opt$out, '.chr',i,'.shapeit_log.snp.mm'))
	system(paste0('rm ',opt$out, '.chr',i,'.shapeit_log.ind.mm'))
	system(paste0('rm ',opt$out, '.chr',i,'.shapeit_log.log'))
	system(paste0('rm ',opt$out,'.chr',i,'.harmonised_temp.haps'))
	system(paste0('rm ',opt$out,'.chr',i,'.harmonised_temp.sample'))
	system(paste0('rm ',opt$out,'.chr',i,'_exclusions'))
	system(paste0('rm ',opt$out,'.chr',i,'.ped'))
	system(paste0('rm ',opt$out,'.chr',i,'.map'))
	system(paste0('rm ',opt$out,'.chr',i,'.log'))

}	

sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
cat('Done!\n')
sink()

###########################################
# Impute genetic data
###########################################

sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
cat('Imputing data using Impute2...')
sink()

system(paste0('mkdir ',opt$out,'.impute2.logs'))

# Look up length of each chromosome
maxPos_table<-foreach(i=CHROMS, .combine=rbind) %dopar% {
	if(i == 'X'){
		data.frame(CHR=i,maxPos=system(paste0('zcat ',opt$ref,'/1000GP_Phase3_chr',i,"_NONPAR.legend.gz | tail -n 1 | cut -d ' ' -f 2"),intern=T), stringsAsFactors=F)
	} else {
		data.frame(CHR=i,maxPos=system(paste0('zcat ',opt$ref,'/1000GP_Phase3_chr',i,".legend.gz | tail -n 1 | cut -d ' ' -f 2"),intern=T), stringsAsFactors=F)
	}
}

# Create a data.frame containing list all 5Mb chunks for imputation
chunks<-NULL
for(i in CHROMS){
  starts<-seq(0,as.numeric(maxPos_table$maxPos[maxPos_table$CHR==i]),5e6)
  ends <- starts+5e6
  chunks<-rbind(chunks,data.frame(CHR=i, Start=starts, End=ends))
}

foreach(i=1:dim(chunks)[1], .combine=c) %dopar% {
	if(chunks$CHR[i] == 'X'){
    imp_log<-system(paste0(opt$impute2,' -m ',opt$ref,'/genetic_map_chr',chunks$CHR[i],'_nonPAR_combined_b37.txt -h ',opt$ref,'/1000GP_Phase3_chr',chunks$CHR[i],'_NONPAR.hap.gz -l ',opt$ref,'/1000GP_Phase3_chr',chunks$CHR[i],'_NONPAR.legend.gz -known_haps_g ',opt$out,'.chr',chunks$CHR[i],'.harmonised.haps -int ',chunks$Start[i],' ',chunks$End[i],' -Ne 20000 -o ',opt$out,'.chr',chunks$CHR[i],'.',chunks$Start[i],'_',chunks$End[i]))
	} else {
		imp_log<-system(paste0(opt$impute2,' -m ',opt$ref,'/genetic_map_chr',chunks$CHR[i],'_combined_b37.txt -h ',opt$ref,'/1000GP_Phase3_chr',chunks$CHR[i],'.hap.gz -l ',opt$ref,'/1000GP_Phase3_chr',chunks$CHR[i],'.legend.gz -known_haps_g ',opt$out,'.chr',chunks$CHR[i],'.harmonised.haps -int ',chunks$Start[i],' ',chunks$End[i],' -Ne 20000 -o ',opt$out,'.chr',chunks$CHR[i],'.',chunks$Start[i],'_',chunks$End[i]))
	}
    #test for memory-lack bug (step_7_log will be 137 if killed, otherwise 0)
    if(imp_log == 137){
			sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
			cat("Chr",chunks$CHR[i],chunks$Start[i],'-',chunks$End[i],":There was a memory problem.\n")
			sink()
	}
}

foreach(i=CHROMS, .combine=c) %dopar% {
  # Combine the chunks into per chromosome files
  starts<-chunks$Start[chunks$CHR == i]
  ends<-chunks$End[chunks$CHR == i]
  chunk_files<-paste0(opt$out,'.chr',i,'.',starts[1:length(starts)],'_',ends[1:length(starts)])
  # Check which chunk_files exist
  chunk_files_present<-NULL
  for(k in 1:length(chunk_files)){
		if(file.exists(chunk_files[k])){
			chunk_files_present<-c(chunk_files_present,chunk_files[k])
		}
  }
  
  system(paste0('cat ',paste(chunk_files_present,collapse=' '),' > ',opt$out,'.chr',i,'.gen'))
  
	# Delete unimportant per chromosome files, and move log files to a log file folder
  for(k in 1:length(chunk_files)){
		system(paste0('rm ',chunk_files[k],'_diplotype_ordering'))
		system(paste0('rm ',chunk_files[k],'_info'))
		system(paste0('rm ',chunk_files[k],'_info_by_sample'))
		system(paste0('rm ',chunk_files[k]))
		system(paste0('mv ',chunk_files[k],'_summary ',opt$out,'.impute2.logs/'))
		system(paste0('mv ',chunk_files[k],'_warnings ',opt$out,'.impute2.logs/'))
  }
  
	system(paste0('rm ',opt$out,'.chr',i,'.harmonised.haps'))
	system(paste0('rm ',opt$out,'.chr',i,'.harmonised.sample'))
}

sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
cat("Done!\n")
sink()

sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
cat("Reformatting...")
sink()

# Update rsids, and save as a plink file and remove missing SNPs
foreach(i=CHROMS, .combine=c) %dopar% {
	# Insert chromosome number
	gen<-fread(paste0(opt$out,'.chr',i,'.gen'))
	gen$V1<-as.character(i)
	
	# Exclude INDELS
	gen<-gen[nchar(gen$V4) == 1 & nchar(gen$V5) == 1,]
	
	# Include only the RSID if available
	gen$V2[grepl('rs',gen$V2)]<-gsub(':.*','',gen$V2[grepl('rs',gen$V2)])
	fwrite(gen, paste0(opt$out,'.chr',i,'.gen'), col.names=F, row.names=F, quote=F, sep=' ')
	rm(gen)

	# Save as a plink file
	system(paste0(opt$qctool,' -g ',opt$out,'.chr',i,'.gen -threshold ',opt$hardcall_thresh,' -og ', opt$out,'.chr',i,' -ofiletype binary_ped'))
	system(paste0('cp ',opt$out,'.fam ', opt$out,'.chr',i,'.fam'))

	# Remove missing SNPs in the plink files
	system(paste0(opt$plink,' --bfile ',opt$out,'.chr',i,' --geno 0.1 --make-bed --out ',opt$out,'.1KGphase3.chr',i,' --memory 2000'))
	
	# Insert the chromosome number to the bim file
	bim<-fread(paste0(opt$out,'.1KGphase3.chr',i,'.bim'))
	bim$V1<-i
	fwrite(bim, paste0(opt$out,'.1KGphase3.chr',i,'.bim'), col.names=F, row.names=F, quote=F, sep='\t')
		
	# Rename the gen file
	system(paste0('mv ',opt$out,'.chr',i,'.gen ',opt$out,'.1KGphase3.chr',i,'.gen'))
	system(paste0('rm ',opt$out,'.chr',i,'.fam'))
	system(paste0('rm ',opt$out,'.chr',i,'.bim'))
	system(paste0('rm ',opt$out,'.chr',i,'.bed'))

}

# Update FID and IID in the fam files
for(i in CHROMS){
		fam<-data.frame(fread(paste0(opt$out,'.1KGphase3.chr',i,'.fam')))
		fam$V1<-opt$FID
		fam$V2<-opt$IID
		write.table(fam,paste0(opt$out,'.1KGphase3.chr',i,'.fam'), col.names=F, row.names=F, quote=F)
}

# Create a sample file for the bgen files
sample_file<-data.frame(ID=c(0,paste0(opt$FID,'_',opt$IID)),
						sex=c('D','M'))			

write.table(sample_file, paste0(opt$out,'.sample'), col.names=T, row.names=F, quote=F)

sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
cat("Done!\n")
sink()

# Delete the log file folders
system(paste0('rm -r ',opt$out,'.impute2.logs'))

# Delete the remaining temporary files
system(paste0('rm ',opt$out,'.ped'))
system(paste0('rm ',opt$out,'.map'))
system(paste0('rm ',opt$out,'.log'))
system(paste0('rm ',opt$out,'.hh'))
system(paste0('rm ',opt$out,'_duplicates.snplist'))

# Count the number of variants after imputation in gen files and plink files
final_n_gen<-system(paste0('wc -l ',opt$out,'*.1KGphase3.chr*.gen'), intern=T)
final_n_gen<-final_n_gen[length(final_n_gen)]
final_n_gen<-gsub(' ','',gsub('total','', final_n_gen))

final_n_plink<-system(paste0('wc -l ',opt$out,'*.1KGphase3.chr*.bim'), intern=T)
final_n_plink<-final_n_plink[length(final_n_plink)]
final_n_plink<-gsub(' ','',gsub('total','', final_n_plink))

sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
cat(paste0('.gen files contain ',final_n_gen,' SNPs.\n'))
cat(paste0('.plink files contain ',final_n_plink,' SNPs.\n'))
sink()

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$out,'.23andMe_imputer.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
