#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()

source('/users/k1806347/brc_scratch/Software/MyGit/GenoPred/config_used/Target_scoring.config')

# Create output directory if it doesn't exists already
system(paste0('mkdir -p ',TEDS_output_dir,'/Genotype/Harmonised'))

library(data.table)

###################
# Identify list of SNPs in the target sample that are in the reference
###################

# Read in target SNP data, and update bim file to have rsIDs
target<-NULL
for(i in 1:22){
	target<-rbind(target, fread(paste0(TEDS_original_dir,'/TEDS.hrc.oee.QC.affy.QC.inclNonOverlapping.chr',i,'.bim')))
	print(i)
}
target<-target[,c('V1','V2','V4','V5','V6')]
names(target)<-c('CHR','SNP','BP','A1','A2')
target$ID<-target$SNP
target$ID[grepl('rs', target$SNP)]<-gsub(':.*','',target$SNP[grepl('rs', target$SNP)])

# Read in reference SNP data
ref<-NULL
for(i in 1:22){
	ref<-rbind(ref, fread(paste0(Geno_1KG_dir,'/1KGPhase3.w_hm3.chr',i,'.bim')))
	print(i)
}
ref<-ref[,c('V1','V2','V4','V5','V6')]
names(ref)<-c('CHR','SNP','BP','A1','A2')

# Merge reference and target by SNP id
ref_target<-merge(ref,target,by.x='SNP', by.y='ID')

# Write a list of SNPs to be extracted from the target sample
fwrite(ref_target[,'SNP.y'], paste0(TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.snplist'), col.names=F)

rm(ref,target,ref_target)
gc()

###################
# Extract SNPs in reference and convert to PLINK format
###################

# Extract SNPs in reference
for(i in 1:22){
	system(paste0(plink1_9,' --bfile ',TEDS_original_dir,'/TEDS.hrc.oee.QC.affy.QC.inclNonOverlapping.chr',i,' --extract ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.snplist --make-bed --out ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.chr',i,' --memory 7000'))
}

##################
# Harmonise alleles with reference
##################

target<-NULL
for(i in 1:22){
target<-rbind(target, fread(paste0(TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.chr',i,'.bim')))
print(i)
}
names(target)<-c('CHR','SNP','POS','BP','A1','A2')
target$ID<-target$SNP
target$ID[grepl('rs', target$SNP)]<-gsub(':.*','',target$SNP[grepl('rs', target$SNP)])

# Read in reference data
ref<-NULL
for(i in 1:22){
ref<-rbind(ref, fread(paste0(Geno_1KG_dir,'/1KGPhase3.w_hm3.chr',i,'.bim')))
print(i)
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
target$SNP_IUPAC<-paste0(target$ID,':',target$IUPAC)

# Create IUPAC codes in ref data
ref$IUPAC[ref$A1 == 'A' & ref$A2 =='T' | ref$A1 == 'T' & ref$A2 =='A']<-'W'
ref$IUPAC[ref$A1 == 'C' & ref$A2 =='G' | ref$A1 == 'G' & ref$A2 =='C']<-'S'
ref$IUPAC[ref$A1 == 'A' & ref$A2 =='G' | ref$A1 == 'G' & ref$A2 =='A']<-'R'
ref$IUPAC[ref$A1 == 'C' & ref$A2 =='T' | ref$A1 == 'T' & ref$A2 =='C']<-'Y'
ref$IUPAC[ref$A1 == 'G' & ref$A2 =='T' | ref$A1 == 'T' & ref$A2 =='G']<-'K'
ref$IUPAC[ref$A1 == 'A' & ref$A2 =='C' | ref$A1 == 'C' & ref$A2 =='A']<-'M'
ref$SNP_IUPAC<-paste0(ref$SNP,':',ref$IUPAC)

# Merge target and reference based on SNP id
ref_target<-merge(ref, target, by.x='SNP', by.y='ID')

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

sink(file = paste0(TEDS_output_dir,'/Genotype/Harmonised/TEDS.harmonisation.log'), append = F)
cat(dim(flip)[1],' variants need to be flipped.\n',sep='')
sink()

target_in_ref_n<-dim(incl)[1]
ref_n<-dim(ref)[1]

sink(file = paste0(TEDS_output_dir,'/Genotype/Harmonised/TEDS.harmonisation.log'), append = T)
cat(target_in_ref_n/ref_n*100,'% of variants in reference are in the target.\n',sep='')
sink()

target$SNP[(!target$SNP_IUPAC %in% incl$SNP_IUPAC.y)]<-paste0(target$SNP[(!target$SNP_IUPAC %in% incl$SNP_IUPAC.y)], '_excl')

target$SNP_IUPAC<-NULL
target$IUPAC<-NULL
target$ID<-NULL

fwrite(target[!grepl('_excl', target$SNP),][,'SNP'], paste0(TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.extract'), col.names=F)
for(i in 1:22){
	fwrite(target[target$CHR == i,], paste0(TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.chr',i,'.bim'), col.names=F, sep=' ')
}

if(dim(flip)[1] > 0){
  fwrite(flip[,'SNP'], paste0(TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.flip'), col.names=F)
}

rm(target,ref,ref_target,incl,flip,flip_tmp)
gc()

for(i in 1:22){
	system(paste0(plink1_9,' --bfile ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.chr',i,' --extract ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.extract --make-bed --out ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.QCd.chr',i,' --memory 7000'))
}

for(i in 1:22){
system(paste0('rm ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.chr',i,'.bed'))
system(paste0('rm ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.chr',i,'.bim'))
system(paste0('rm ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.chr',i,'.fam'))
system(paste0('rm ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.chr',i,'.log'))
}                                    
system(paste0('rm ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.extract'))
system(paste0('rm ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.snplist'))

# Change SNP IDs to contain only the rsID
for(i in 1:22){
target<-fread(paste0(TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.QCd.chr',i,'.bim'))
target$V2[grepl('rs', target$V2)]<-gsub(':.*','',target$V2[grepl('rs', target$V2)])
write.table(target, paste0(TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.QCd.chr',i,'.bim'), col.names=F, row.names=F, quote=F)
}

##################
# Insert missing SNPs into the reference data
##################

sink(file = paste0(TEDS_output_dir,'/Genotype/Harmonised/TEDS.harmonisation.log'), append = T)
cat('Inserting missing HapMap3 SNPs...',sep='')
sink()

# Update IDs in reference to avoid conflict with the target
ref_fam<-fread(paste0(Geno_1KG_dir,'/1KGPhase3.w_hm3.chr22.fam'))
ref_ID_update<-data.frame(ref_fam$V1, ref_fam$V2, paste0(ref_fam$V1,'_REF'),paste0(ref_fam$V2,'_REF'))
fwrite(ref_ID_update, paste0(TEDS_output_dir,'/Genotype/Harmonised/ref_ID_update.txt'), sep=' ', col.names=F)

for(chr in 1:22){
  system(paste0(plink1_9,' --bfile ',Geno_1KG_dir,'/1KGPhase3.w_hm3.chr',chr,' --make-bed --update-ids ',TEDS_output_dir,'/Genotype/Harmonised/ref_ID_update.txt --out ',TEDS_output_dir,'/Genotype/Harmonised/REF.chr',chr,' --memory 7000'))

  # Merge target and reference plink files to insert missing SNPs
  system(paste0(plink1_9,' --bfile ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.QCd.chr',chr,' --bmerge ',TEDS_output_dir,'/Genotype/Harmonised/REF.chr',chr,' --make-bed --out ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.QCd.tmp.chr',chr,' --memory 7000'))

  # Extract only target individuals
  system(paste0(plink1_9,' --bfile ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.QCd.tmp.chr',chr,' --remove ',TEDS_output_dir,'/Genotype/Harmonised/REF.chr',chr,'.fam --make-bed --out ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.QCd.AllSNP.chr',chr,' --memory 7000'))
  system(paste0('rm ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.QCd.tmp.chr',chr,'*'))
  system(paste0('rm ',TEDS_output_dir,'/Genotype/Harmonised/REF.chr',chr,'*'))
  system(paste0('rm ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.QCd.AllSNP.chr',chr,'.nosex'))
  system(paste0('rm ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.QCd.AllSNP.chr',chr,'.log'))
}
system(paste0('rm ',TEDS_output_dir,'/Genotype/Harmonised/ref_ID_update.txt'))

system(paste0('rm ',TEDS_output_dir,'/Genotype/Harmonised/TEDS.w_hm3.QCd.chr*'))

sink(file = paste0(TEDS_output_dir,'/Genotype/Harmonised/TEDS.harmonisation.log'), append = T)
cat('Done!\n')
sink()
