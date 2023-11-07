#!/usr/bin/Rscript

# To speed up the basic implementation of the pipeline, I am going to store a pre-processed version of 1KG on zenodo
# The 1KG data has been prepared using the previously used prep_1kg package.
# We are going to change the layout of the reference data, so it is more easily replicated using other reference datasets

# <ref_dir>
# ├── ref.chr<1-22>.<bed/bim/fam>
# ├── ref.pop.txt # FID, IID and POP (with header)
# ├── ref.keep.list # POP and PATH (without header)
# ├── keep_files
# │   └──<POP>.keep # FID and IID (without header)
# └── freq_files
#     └──<POP>
#         └──ref.<POP>.chr<CHR>.frq # PLINK .frq format

# Create a folder to store the preprocessed 1kg data in
setwd('~/oliverpainfel/Software/MyGit/GenoPred/GenoPredPipe/')
dir.create('1kg')

# Copy plink bed bim fam files over
for( i in 1:22 ){
    for( j in c('bed','bim','fam') ){
        system(paste0('cp resources/data/1kg/1KGPhase3.w_hm3.chr',i,'.',j,' 1kg/ref.chr',i,'.',j))
    }
}

# Format and cp population data
library(data.table)
pop_dat<-fread('resources/data/1kg/ref_pop_dat.txt')
pop_dat<-data.frame(
    FID=pop_dat$sample,
    IID=pop_dat$sample,
    POP=pop_dat$super_pop
)
write.table(pop_dat, '1kg/ref.pop.txt', col.names=T, row.names=F, quote=F)

# Copy over keep and freq files
dir.create('1kg/keep_files')
for(i in unique(pop_dat$POP)){
    system(paste0('cp resources/data/1kg/keep_files/',i,'_samples.keep 1kg/keep_files/',i,'.keep'))
}

dir.create('1kg/freq_files')
for(i in unique(pop_dat$POP)){
    dir.create(paste0('1kg/freq_files/',i))
    for(j in 1:22){
        system(paste0('cp resources/data/1kg/freq_files/',i,'/1KGPhase3.w_hm3.',i,'.chr',j,'.frq 1kg/freq_files/', i,'/ref.',i,'.chr',j,'.frq'))
    }
}

# Copy over super_pop_keep.list
keep_list<-fread('resources/data/1kg/super_pop_keep.list', header=F)
keep_list$V2<-gsub('1kg','ref',keep_list$V2)
keep_list$V2<-gsub('_samples','',keep_list$V2)
write.table(keep_list, '1kg/ref.keep.list', col.names=F, row.names=F, quote=F)

# Compress folder and upload to zenodo
system(paste0('tar -czvf genopredpipe_1kg.tar.gz 1kg'))
