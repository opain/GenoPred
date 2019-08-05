#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
	make_option("--chr", action="store", default=NA, type='numeric',
			help="Chromosome number [required]"),
	make_option("--input_dir", action="store", default=NA, type='character',
			help="Directory containing imputed UKBB data [required]"),
	make_option("--target_fam", action="store", default=NA, type='character',
			help="Path to fam file for target sample [required]"),
	make_option("--reference_dir", action="store", default=NA, type='character',
			help="Directory containing reference data [required]"),
	make_option("--qctool2", action="store", default=NA, type='character',
			help="Path to QCTOOL V2 binary [required]"),
	make_option("--plink", action="store", default=NA, type='character',
			help="Path to PLINK V1.9 binary [required]"),		
	make_option("--output_dir", action="store", default=NA, type='character',
			help="Output folder name [required]"),
	make_option("--debug", action="store", default=F, type='logical',
			help="Set to T to create object for debugging [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Set exports to avoid core hungry software
system(paste0('export MKL_NUM_THREADS=1'))
system(paste0('export NUMEXPR_NUM_THREADS=1'))
system(paste0('export OMP_NUM_THREADS=1'))

# Create output directory if it doesn't exists already
system(paste0('mkdir -p ',opt$output_dir))

if(opt$debug == T){
	save(opt,file=paste0(opt$output_dir,'/Harmonisation_of_UKBB.debug.opt'))
	q()
}

library(data.table)

chr<-opt$chr

###################
# Identify list of SNPs in the target sample that are in the reference
###################

# Read in target SNP data
target<-fread(paste0(opt$input_dir,'/ukb_mfi_chr',chr,'_v3.txt'))
target<-target[,c('V2','V3','V4','V5')]
names(target)<-c('SNP','BP','A1','A2')

# Read in reference SNP data
ref<-fread(paste0(opt$reference_dir,'/1KGPhase3.w_hm3.chr',chr,'.bim'))
ref<-ref[,c('V2','V4','V5','V6')]
names(ref)<-c('SNP','BP','A1','A2')

# Merge reference and target by SNP id
ref_target<-merge(ref,target,by='SNP')

# Write a list of SNPs to be extracted from the target sample
fwrite(ref_target[,'SNP'], paste0(opt$output_dir,'/UKBB.w_hm3.chr',chr,'.snplist'), col.names=F)

rm(target,ref,ref_target)
gc()

###################
# Extract SNPs in reference and convert to PLINK format
###################

# Extract SNPs in reference and convert to plink using qctool2
system(paste0(opt$qctool2,' -g ',opt$input_dir,'/ukb_imp_chr',chr,'_v3.bgen -s ',opt$target_fam,' -incl-rsids ',opt$output_dir,'/UKBB.w_hm3.chr',chr,'.snplist -ofiletype binary_ped -og ',opt$output_dir,'/UKBB.w_hm3.chr',chr,''))

##################
# Harmonise alleles with reference
##################

# Read in target data
target<-fread(paste0(opt$output_dir,'/UKBB.w_hm3.chr',chr,'.bim'))
names(target)<-c('CHR','SNP','POS','BP','A1','A2')

# Read in reference data
ref<-fread(paste0(opt$reference_dir,'/1KGPhase3.w_hm3.chr',chr,'.bim'))
ref<-ref[,c('V2','V4','V5','V6')]
names(ref)<-c('SNP','BP','A1','A2')

# Create IUPAC codes in target data
target$IUPAC[target$A1 == 'A' & target$A2 =='T' | target$A1 == 'T' & target$A2 =='A']<-'W'
target$IUPAC[target$A1 == 'C' & target$A2 =='G' | target$A1 == 'G' & target$A2 =='C']<-'S'
target$IUPAC[target$A1 == 'A' & target$A2 =='G' | target$A1 == 'G' & target$A2 =='A']<-'R'
target$IUPAC[target$A1 == 'C' & target$A2 =='T' | target$A1 == 'T' & target$A2 =='C']<-'Y'
target$IUPAC[target$A1 == 'G' & target$A2 =='T' | target$A1 == 'T' & target$A2 =='G']<-'K'
target$IUPAC[target$A1 == 'A' & target$A2 =='C' | target$A1 == 'C' & target$A2 =='A']<-'M'
target$SNP_IUPAC<-paste0(target$SNP,':',target$IUPAC)

# Create IUPAC codes in target data
ref$IUPAC[ref$A1 == 'A' & ref$A2 =='T' | ref$A1 == 'T' & ref$A2 =='A']<-'W'
ref$IUPAC[ref$A1 == 'C' & ref$A2 =='G' | ref$A1 == 'G' & ref$A2 =='C']<-'S'
ref$IUPAC[ref$A1 == 'A' & ref$A2 =='G' | ref$A1 == 'G' & ref$A2 =='A']<-'R'
ref$IUPAC[ref$A1 == 'C' & ref$A2 =='T' | ref$A1 == 'T' & ref$A2 =='C']<-'Y'
ref$IUPAC[ref$A1 == 'G' & ref$A2 =='T' | ref$A1 == 'T' & ref$A2 =='G']<-'K'
ref$IUPAC[ref$A1 == 'A' & ref$A2 =='C' | ref$A1 == 'C' & ref$A2 =='A']<-'M'
ref$SNP_IUPAC<-paste0(ref$SNP,':',ref$IUPAC)

# Merge target and reference based on SNP id
ref_target<-merge(ref, target, by='SNP')

# Identify SNPs for which alleles need to be flipped
flip_tmp<-ref_target[(ref_target$IUPAC.x == 'R' & ref_target$IUPAC.y == 'Y' | 
                      ref_target$IUPAC.x == 'Y' & ref_target$IUPAC.y == 'R' | 
                      ref_target$IUPAC.x == 'K' & ref_target$IUPAC.y == 'M' |
                      ref_target$IUPAC.x == 'M' & ref_target$IUPAC.y == 'K'),]

# Idenitfy SNPs which match the reference alleles
incl<-ref_target[ ref_target$IUPAC.x == 'R' & ref_target$IUPAC.y == 'R' | 
                  ref_target$IUPAC.x == 'Y' & ref_target$IUPAC.y == 'Y' | 
                  ref_target$IUPAC.x == 'K' & ref_target$IUPAC.y == 'K' |
                  ref_target$IUPAC.x == 'M' & ref_target$IUPAC.y == 'M' ]

# If a SNP that needs to be flipped has a duplicate that is on the correct strand, remove it.
flip<-flip_tmp[!(flip_tmp$SNP %in% incl$SNP)]

# Combine SNPs that match and those that need to be flipped.
incl<-rbind(incl,flip)

sink(file = paste(opt$output_dir,'/UKBB.harmonisation.chr',chr,'.log',sep=''), append = T)
cat(dim(flip)[1],' variants need to be flipped on chromosome ',chr,'.\n',sep='')
sink()

target_in_ref_n<-dim(incl)[1]
ref_n<-dim(ref)[1]

sink(file = paste(opt$output_dir,'/UKBB.harmonisation.chr',chr,'.log',sep=''), append = T)
cat(target_in_ref_n/ref_n*100,'% variants on ',chr,' in reference are in the target.\n',sep='')
sink()

target$SNP[(!target$SNP_IUPAC %in% incl$SNP_IUPAC.y)]<-paste0(target$SNP[(!target$SNP_IUPAC %in% incl$SNP_IUPAC.y)], '_excl')

target$SNP_IUPAC<-NULL
target$IUPAC<-NULL

fwrite(target[!grepl('_excl', target$SNP),][,'SNP'], paste0(opt$output_dir,'/UKBB.w_hm3.chr',chr,'.extract'), col.names=F)
fwrite(target, paste0(opt$output_dir,'/UKBB.w_hm3.chr',chr,'.bim'), col.names=F, sep=' ')

if(dim(flip)[1] > 0){
  fwrite(flip[,'SNP'], paste0(opt$output_dir,'/UKBB.w_hm3.chr',chr,'.flip'), col.names=F)
}

rm(target,ref,ref_target,incl,flip,flip_tmp)
gc()

system(paste0(opt$plink,' --bfile ',opt$output_dir,'/UKBB.w_hm3.chr',chr,' --extract ',opt$output_dir,'/UKBB.w_hm3.chr',chr,'.extract --make-bed --out ',opt$output_dir,'/UKBB.w_hm3.QCd.chr',chr,' --memory 7000'))

system(paste0('rm ',opt$output_dir,'/UKBB.w_hm3.chr',chr,'.*'))
system(paste0('rm ',opt$output_dir,'/UKBB.w_hm3.QCd.chr',chr,'.nosex'))

##################
# Insert missing SNPs into the reference data
##################

sink(file = paste(opt$output_dir,'/UKBB.harmonisation.chr',chr,'.log',sep=''), append = T)
cat('Inserting missing HapMap3 SNPs...',sep='')
sink()

# Update IDs in reference to avoid conflict with the target
ref_fam<-fread(paste0(opt$reference_dir,'/1KGPhase3.w_hm3.chr22.fam'))
ref_ID_update<-data.frame(ref_fam$V1, ref_fam$V2, paste0(ref_fam$V1,'_REF'),paste0(ref_fam$V2,'_REF'))
fwrite(ref_ID_update, paste0(opt$output_dir,'/ref_ID_update_chr',chr,'.txt'), sep=' ', col.names=F)
system(paste0(opt$plink,' --bfile ',opt$reference_dir,'/1KGPhase3.w_hm3.chr',chr,' --make-bed --update-ids ',opt$output_dir,'/ref_ID_update_chr',chr,'.txt --out ',opt$output_dir,'/REF.chr',chr,' --memory 7000'))

# Merge target and reference plink files to insert missing SNPs
system(paste0(opt$plink,' --bfile ',opt$output_dir,'/UKBB.w_hm3.QCd.chr',chr,' --bmerge ',opt$output_dir,'/REF.chr',chr,' --make-bed --out ',opt$output_dir,'/UKBB.w_hm3.QCd.tmp.chr',chr,' --memory 7000'))
system(paste0('rm ',opt$output_dir,'/UKBB.w_hm3.QCd.chr',chr,'*'))

# Extract only target individuals
system(paste0(opt$plink,' --bfile ',opt$output_dir,'/UKBB.w_hm3.QCd.tmp.chr',chr,' --remove ',opt$output_dir,'/REF.chr',chr,'.fam --make-bed --out ',opt$output_dir,'/UKBB.w_hm3.QCd.AllSNP.chr',chr,' --memory 7000'))
system(paste0('rm ',opt$output_dir,'/UKBB.w_hm3.QCd.tmp.chr',chr,'*'))
system(paste0('rm ',opt$output_dir,'/REF.chr',chr,'*'))
system(paste0('rm ',opt$output_dir,'/UKBB.w_hm3.QCd.AllSNP.chr',chr,'.nosex'))
system(paste0('rm ',opt$output_dir,'/UKBB.w_hm3.QCd.AllSNP.chr',chr,'.log'))
system(paste0('rm ',opt$output_dir,'/ref_ID_update_chr',chr,'.txt'))

sink(file = paste(opt$output_dir,'/UKBB.harmonisation.chr',chr,'.log',sep=''), append = T)
cat('Done!\n')
sink()






