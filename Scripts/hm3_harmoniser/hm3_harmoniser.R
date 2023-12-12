#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--target", action="store", default=NA, type='character',
		help="Path to per chromosome target sample plink files [required]"),
make_option("--ref", action="store", default=NA, type='character',
		help="Path to per chromosome target sample plink files [required]"),
make_option("--plink", action="store", default=NA, type='character',
		help="Path to plink1.9 [required]"),
make_option("--out", action="store", default=NA, type='character',
		help="Path for output files [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

opt$out_dir<-paste0(dirname(opt$out),'/')
system(paste0('mkdir -p ',opt$out_dir))

sink(file = paste(opt$out,'.hm3_harmoniser.log',sep=''), append = F)
cat(
"#################################################################
# hm3_harmoniser.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at",as.character(start.time),'
Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

library(data.table)

###################
# Identify list of SNPs in the target sample that are in the reference
###################

sink(file = paste(opt$out,'.hm3_harmoniser.log',sep=''), append = T)
cat('Extracting HapMap3 SNPs...')
sink()

# Read in reference SNP data
ref<-NULL
for(i in 1:22){
  ref<-rbind(ref, fread(paste0(opt$ref,i,'.bim')))
}
ref<-ref[,c('V1','V2','V4','V5','V6')]
names(ref)<-c('CHR','SNP','BP','A1','A2')

# Read in target SNP data
target<-NULL
for(i in 1:22){
  target_tmp<-fread(paste0(opt$target, i,'.bim'))
  target_tmp<-target_tmp[,c('V1','V2','V4','V5','V6')]
  names(target_tmp)<-c('CHR','SNP','BP','A1','A2')
  target_tmp<-target_tmp[target_tmp$SNP %in% ref$SNP]
  target<-rbind(target, target_tmp)
}

# Write a list of SNPs to be extracted from the target sample
fwrite(target[,'SNP'], paste0(opt$out,'.w_hm3.snplist'), col.names=F)

rm(target,ref)
gc()

###################
# Extract SNPs in reference and convert to PLINK format
###################

# Extract SNPs in reference
for(i in 1:22){
  system(paste0(opt$plink,' --bfile ',opt$target,i,' --extract ',opt$out,'.w_hm3.snplist --make-bed --out ',opt$out,i))
}

sink(file = paste(opt$out,'.hm3_harmoniser.log',sep=''), append = T)
cat('Done!\n')
sink()

##################
# Harmonise alleles with reference
##################

sink(file = paste(opt$out,'.hm3_harmoniser.log',sep=''), append = T)
cat('Harmonising alleles with reference...\n')
sink()

# Read in target data
target<-NULL
for(i in 1:22){
  target<-rbind(target, fread(paste0(opt$out,i,'.bim')))
}
target<-target[,c('V1','V2','V4','V5','V6')]
names(target)<-c('CHR','SNP','BP','A1','A2')

# Read in reference data
ref<-NULL
for(i in 1:22){
  ref<-rbind(ref, fread(paste0(opt$ref,i,'.bim')))
}
ref<-ref[,c('V1','V2','V4','V5','V6')]
names(ref)<-c('CHR','SNP','BP','A1','A2')

# Create IUPAC codes in target data
target$IUPAC[target$A1 == 'A' & target$A2 =='T' | target$A1 == 'T' & target$A2 =='A']<-'W'
target$IUPAC[target$A1 == 'C' & target$A2 =='G' | target$A1 == 'G' & target$A2 =='C']<-'S'
target$IUPAC[target$A1 == 'A' & target$A2 =='G' | target$A1 == 'G' & target$A2 =='A']<-'R'
target$IUPAC[target$A1 == 'C' & target$A2 =='T' | target$A1 == 'T' & target$A2 =='C']<-'Y'
target$IUPAC[target$A1 == 'G' & target$A2 =='T' | target$A1 == 'T' & target$A2 =='G']<-'K'
target$IUPAC[target$A1 == 'A' & target$A2 =='C' | target$A1 == 'C' & target$A2 =='A']<-'M'
target$SNP_IUPAC<-paste0(target$SNP,':',target$IUPAC)

# Create IUPAC codes in ref data
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

sink(file = paste(opt$out,'.hm3_harmoniser.log',sep=''), append = T)
cat(dim(flip)[1],' variants need to be flipped.\n',sep='')
sink()

target_in_ref_n<-nrow(incl)
ref_n<-nrow(ref)

sink(file = paste(opt$out,'.hm3_harmoniser.log',sep=''), append = T)
cat(target_in_ref_n," of ", ref_n," reference variants (", round(target_in_ref_n/ref_n*100,2), "%) are in the target.\n",sep='')
sink()

target$SNP[(!target$SNP_IUPAC %in% incl$SNP_IUPAC.y)]<-paste0(target$SNP[(!target$SNP_IUPAC %in% incl$SNP_IUPAC.y)], '_excl')

target$SNP_IUPAC<-NULL
target$IUPAC<-NULL

fwrite(target[!grepl('_excl', target$SNP),][,'SNP'], paste0(opt$out,'.w_hm3.extract'), col.names=F)

target$POS<-0
target<-target[,c('CHR','SNP','POS','BP','A1','A2')]

for(i in 1:22){
  fwrite(target[target$CHR == i,], paste0(opt$out,i,'.bim'), col.names=F, sep=' ')
}

if(dim(flip)[1] > 0){
  fwrite(flip[,'SNP'], paste0(opt$out,'.flip'), col.names=F)
}

rm(target,ref,ref_target,incl,flip,flip_tmp)
gc()

for(i in 1:22){
  system(paste0(opt$plink,' --bfile ',opt$out,i,' --extract ',opt$out,'.w_hm3.extract --make-bed --out ',opt$out,i))
}

system(paste0('rm ',opt$out,'*~'))
system(paste0('rm ',opt$out,'.snplist'))
system(paste0('rm ',opt$out,'extract'))

##################
# Insert missing SNPs into the reference data
##################

sink(file = paste(opt$out,'.hm3_harmoniser.log',sep=''), append = T)
cat('Inserting missing HapMap3 SNPs...',sep='')
sink()

# Update IDs in reference to avoid conflict with the target
ref_fam<-fread(paste0(opt$ref,'22.fam'))
ref_ID_update<-data.frame(ref_fam$V1, ref_fam$V2, paste0(ref_fam$V1,'_REF'),paste0(ref_fam$V2,'_REF'))
fwrite(ref_ID_update, paste0(opt$out,'.ref_ID_update.txt'), sep=' ', col.names=F)

for(i in 1:22){
  system(paste0(opt$plink,' --bfile ',opt$ref,i,' --make-bed --update-ids ',opt$out,'.ref_ID_update.txt --out ',opt$out,'.REF.',i,' --memory 7000'))

  # Merge target and reference plink files to insert missing SNPs
  system(paste0(opt$plink,' --bfile ',opt$out,i,' --bmerge ',opt$out,'.REF.',i,' --make-bed --out ',opt$out,'.tmp',i))

  # Extract only target individuals
  system(paste0(opt$plink,' --bfile ',opt$out,'.tmp',i,' --remove ',opt$out,'.REF.',i,'.fam --make-bed --out ',opt$out,i))

  system(paste0('rm ',opt$out,'.tmp',i,'*'))
  system(paste0('rm ',opt$out,'.REF.',i,'*'))
  system(paste0('rm ',opt$out,i,'.nosex'))
  system(paste0('rm ',opt$out,i,'.log'))
}

system(paste0('rm ',opt$out,'.ref_ID_update.txt'))

sink(file = paste(opt$out,'.hm3_harmoniser.log',sep=''), append = T)
cat('Done!\n')
sink()

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$out,'.hm3_harmoniser.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
