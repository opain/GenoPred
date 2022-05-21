#!/usr/bin/Rscript

###################
# PLINK format 1000 genomes phase 3 data
###################
# Based on instructions from Hannah Meyer https://cran.r-project.org/web/packages/plinkQC/vignettes/Genomes1000.pdf
# With adaptations from Joni Coleman's reference resource

output_dir<-'resources/data/1kg'
dir.create(output_dir)

# Download the 1000Genomes data provided by PLINK (PLINK 2 format)
system(paste0('wget -q -o /dev/null -O ',output_dir,'/all_phase3.pgen.zst https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1'))
system(paste0('wget -q -o /dev/null -O ',output_dir,'/all_phase3.pvar.zst https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1'))
system(paste0('wget -q -o /dev/null -O ',output_dir,'/all_phase3.psam https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam?dl=1'))

# Decompress the pgen file
system(paste0('plink2 --zst-decompress ',output_dir,'/all_phase3.pgen.zst > ',output_dir,'/all_phase3.pgen'))

# Convert to plink 1 format
system(paste0('plink2 --pfile ',output_dir,'/all_phase3 vzs --max-alleles 2 --make-bed --out ',output_dir,'/all_phase3'))

# Delete the plink2 format data
system(paste0('rm ',output_dir,'/all_phase3.pgen*'))
system(paste0('rm ',output_dir,'/all_phase3.pvar*'))

# Make FID == IID
fam<-read.table(paste0(output_dir,'/all_phase3.fam'), header=F)
fam$V1<-fam$V2
write.table(fam, paste0(output_dir,'/all_phase3.fam'), col.names=F, row.names=F, quote=F)

# Split the plink data by chromosome (autosome only)
for(chr in 1:22){
  system(paste0('plink --bfile ',output_dir,'/all_phase3 --chr ',chr,' --allow-extra-chr --make-bed --out ',output_dir,'/all_phase3.chr',chr))
}

# Delete genome wide data
system(paste0('rm ',output_dir,'/all_phase3.bed'))
system(paste0('rm ',output_dir,'/all_phase3.bim'))
system(paste0('rm ',output_dir,'/all_phase3.fam'))
system(paste0('rm ',output_dir,'/all_phase3.log'))

####################
# Harmonise with HapMap3 snplist
####################

library(data.table)

hapmap3_snps<-fread('resources/data/hm3_snplist/w_hm3.snplist')

# Generate IUPAC codes
hapmap3_snps$IUPAC[hapmap3_snps$A1 == 'A' & hapmap3_snps$A2 =='T' | hapmap3_snps$A1 == 'T' & hapmap3_snps$A2 =='A']<-'W'
hapmap3_snps$IUPAC[hapmap3_snps$A1 == 'C' & hapmap3_snps$A2 =='G' | hapmap3_snps$A1 == 'G' & hapmap3_snps$A2 =='C']<-'S'
hapmap3_snps$IUPAC[hapmap3_snps$A1 == 'A' & hapmap3_snps$A2 =='G' | hapmap3_snps$A1 == 'G' & hapmap3_snps$A2 =='A']<-'R'
hapmap3_snps$IUPAC[hapmap3_snps$A1 == 'C' & hapmap3_snps$A2 =='T' | hapmap3_snps$A1 == 'T' & hapmap3_snps$A2 =='C']<-'Y'
hapmap3_snps$IUPAC[hapmap3_snps$A1 == 'G' & hapmap3_snps$A2 =='T' | hapmap3_snps$A1 == 'T' & hapmap3_snps$A2 =='G']<-'K'
hapmap3_snps$IUPAC[hapmap3_snps$A1 == 'A' & hapmap3_snps$A2 =='C' | hapmap3_snps$A1 == 'C' & hapmap3_snps$A2 =='A']<-'M'

hapmap3_snps$SNP_IUPAC<-paste(hapmap3_snps$SNP,hapmap3_snps$IUPAC,sep=':')

# For each chr
for(chr in 1:22){
  bim<-fread(paste0(output_dir,'/all_phase3.chr',chr,'.bim'))
  
  bim$IUPAC[bim$V5 == 'A' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='A']<-'W'
  bim$IUPAC[bim$V5 == 'C' & bim$V6 =='G' | bim$V5 == 'G' & bim$V6 =='C']<-'S'
  bim$IUPAC[bim$V5 == 'A' & bim$V6 =='G' | bim$V5 == 'G' & bim$V6 =='A']<-'R'
  bim$IUPAC[bim$V5 == 'C' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='C']<-'Y'
  bim$IUPAC[bim$V5 == 'G' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='G']<-'K'
  bim$IUPAC[bim$V5 == 'A' & bim$V6 =='C' | bim$V5 == 'C' & bim$V6 =='A']<-'M'
  
  bim$SNP_IUPAC<-paste(bim$V2,bim$IUPAC,sep=':')
  
  bim_hapmap3_snps<-merge(hapmap3_snps,bim,by='SNP_IUPAC')
  sum(duplicated(bim_hapmap3_snps$V2))
  
  bim$V2[!(bim$SNP_IUPAC %in% bim_hapmap3_snps$SNP_IUPAC)] <- paste0(bim$V2[!(bim$SNP_IUPAC %in% bim_hapmap3_snps$SNP_IUPAC)],'_excl')
  bim[!(bim$SNP_IUPAC %in% bim_hapmap3_snps$SNP_IUPAC),]
  
  fwrite(bim[,1:6], paste0(output_dir,'/all_phase3.chr',chr,'.bim'), sep='\t', col.names=F)
  fwrite(as.list(bim_hapmap3_snps$V2), paste0(output_dir,'/all_phase3.chr',chr,'.extract'), sep='\n', col.names=F)   
}

for(chr in 1:22){
  system(paste0('plink --bfile ',output_dir,'/all_phase3.chr',chr,' --extract ',output_dir,'/all_phase3.chr',chr,'.extract --make-bed --out ',output_dir,'/1KGPhase3.w_hm3.chr',chr))
}

####################
# Prepare keep files for super_populations and populations
####################

# Create a keep file listing each population super population from the reference.
dir.create(paste0(output_dir,'/keep_files'))

pop_data<-read.table(paste0(output_dir, '/all_phase3.psam'), header=F, stringsAsFactors=F)

for(i in unique(pop_data$V6)){
  write.table(cbind(pop_data$V1[pop_data$V6 == i],pop_data$V1[pop_data$V6 == i]), paste0(output_dir,'/keep_files/',i,'_samples.keep'), col.names=F, row.names=F, quote=F)
}

for(i in unique(pop_data$V5)){
  write.table(cbind(pop_data$V1[pop_data$V5 == i],pop_data$V1[pop_data$V5 == i]), paste0(output_dir,'/keep_files/',i,'_samples.keep'), col.names=F, row.names=F, quote=F)
}

# Create a file listing the code of each population and the location of the keep file
pop_keep_loc<-data.frame(pop=unique(pop_data$V6))
pop_keep_loc$keep<-paste0(output_dir,'/keep_files/',pop_keep_loc$pop,'_samples.keep')

super_pop_keep_loc<-data.frame(pop=unique(pop_data$V5))
super_pop_keep_loc$keep<-paste0(output_dir, '/keep_files/',super_pop_keep_loc$pop,'_samples.keep')

write.table(super_pop_keep_loc, paste0(output_dir,'/super_pop_keep.list'), col.names=F, row.names=F, quote=F)
write.table(pop_keep_loc, paste0(output_dir,'/pop_keep.list'), col.names=F, row.names=F, quote=F)
write.table(rbind(super_pop_keep_loc,pop_keep_loc), paste0(output_dir,'/super_pop_and_pop_keep.list'), col.names=F, row.names=F, quote=F)

for(i in unique(pop_data$V5)){
  pop_keep_i_loc<-data.frame(pop=unique(pop_data$V6[pop_data$V5 == i]))
  pop_keep_i_loc$keep<-paste0(output_dir,'/keep_files/',pop_keep_i_loc$pop,'_samples.keep')
  write.table(pop_keep_i_loc, paste0(output_dir,'/pop_keep_for_',i,'.list'), col.names=F, row.names=F, quote=F)
}

# Write a file listing ancestry of each reference individual
pop_data<-pop_data[,c('V1','V5','V6')]
names(pop_data)<-c('sample','super_pop','pop')
write.table(pop_data, paste0(output_dir,'/ref_pop_dat.txt'), col.names=T, row.names=F, quote=F)

for(i in unique(pop_data$super_pop)){
  write.table(pop_data[pop_data$super_pop == i,c('sample','pop')], paste0(output_dir,'/ref_pop_dat_for_',i,'.txt'), col.names=T, row.names=F, quote=F)
}

system(paste0('rm ',output_dir,'/all_phase3*'))

####################
# Compute allele frequencies across all individuals
####################

dir.create(paste0(output_dir,'/freq_files'))
for(pop in super_pop_keep_loc$pop){
  dir.create(paste0(output_dir,'/freq_files/',pop))
  for(chr in 1:22){
    system(paste0('plink --bfile ',output_dir,'/1KGPhase3.w_hm3.chr',chr,' --freq --keep ',output_dir,'/keep_files/',pop,'_samples.keep --out ',output_dir,'/freq_files/',pop,'/1KGPhase3.w_hm3.',pop,'.chr',chr))
  }
}

for(pop in pop_keep_loc$pop){
  dir.create(paste0(output_dir,'/freq_files/',pop))
  for(chr in 1:22){
    system(paste0('plink --bfile ',output_dir,'/1KGPhase3.w_hm3.chr',chr,' --freq --keep ',output_dir,'/keep_files/',pop,'_samples.keep --out ',output_dir,'/freq_files/',pop,'/1KGPhase3.w_hm3.',pop,'.chr',chr))
  }
}

dir.create(paste0(output_dir,'/freq_files/AllAncestry'))
system(paste0('plink --bfile ',output_dir,'/1KGPhase3.w_hm3.chr',chr,' --freq --out ',output_dir,'/freq_files/AllAncestry/1KGPhase3.w_hm3.AllAncestry.chr',chr))

# Delete the log and nosex files
system(paste0('rm ',output_dir,'/freq_files/*/1KGPhase3.w_hm3.*.chr*.log'))
