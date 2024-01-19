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

###################
# Prepare summary data for GWAS sumstat cleaning
###################

# Create snp_modifyBuild_offline
make_executable <- function(exe) {
  Sys.chmod(exe, mode = (file.info(exe)$mode | "111"))
}

snp_modifyBuild_offline<-function (info_snp, liftOver, chain, from = "hg18", to = "hg19"){
  tmp_folder<-tempdir()

  if (!all(c("chr", "pos") %in% names(info_snp)))
    stop2("Please use proper names for variables in 'info_snp'. Expected %s.",
          "'chr' and 'pos'")
  liftOver <- normalizePath(liftOver)
  make_executable(liftOver)
  BED <- tempfile(fileext = ".BED")
  info_BED <- with(info_snp, data.frame(paste0("chr", chr),
                                        pos0 = pos - 1L, pos, id = seq_len(nrow(info_snp))))
  fwrite(info_BED, BED, col.names = FALSE, sep = " ")
  lifted<-paste0(tmp_folder,'/tmp.lifted')
  unmapped<-paste0(tmp_folder,'/tmp.unmapped')
  system(paste(liftOver, BED, chain, lifted, unmapped))
  new_pos <- fread(lifted)
  bad <- grep("^#", readLines(unmapped), value = TRUE, invert = TRUE)
  print(paste0(length(bad)," variants have not been mapped."))
  info_snp$pos <- NA
  info_snp$pos[new_pos$V4] <- new_pos$V3
  info_snp
}

output_dir<-'resources/data/ref'
pop_dat<-fread(paste0(output_dir,'/ref.pop.txt'))

chrs<-c(1:22)
for(chr in chrs){

    ######
    # Read in the reference data
    ######
    ref<-list()
    ref[['GRCh37']]<-fread(paste0(output_dir,'/ref.chr',chr,'.bim'))
    ref[['GRCh37']]$V3<-NULL
    names(ref[['GRCh37']])<-c('chr','snp','pos','a1','a2')

    ######
    # Liftover from GRCh37 to GRCh38 and GRCh36
    ######

    # Liftover BP to GRCh38
    ref[['GRCh38']]<-snp_modifyBuild_offline(ref[['GRCh37']], liftOver='/users/k1806347/oliverpainfel/Software/MyGit/GenoDisc/pipeline/resources/software/liftover/liftover', chain='/users/k1806347/oliverpainfel/Software/MyGit/GenoDisc/pipeline/resources/data/liftover/hg19ToHg38.over.chain.gz', from = "hg19", to = "hg38")
    # Liftover BP to GRCh36
    ref[['GRCh36']]<-snp_modifyBuild_offline(ref[['GRCh37']], liftOver='/users/k1806347/oliverpainfel/Software/MyGit/GenoDisc/pipeline/resources/software/liftover/liftover', chain='/users/k1806347/oliverpainfel/Software/MyGit/GenoDisc/pipeline/resources/data/liftover/hg19ToHg18.over.chain.gz', from = "hg19", to = "hg18")

    # Combine the three builds 
    tmp<-ref[['GRCh37']]
    names(tmp)<-c('CHR','SNP','BP_GRCh37','A1','A2')
    tmp$BP_GRCh38<-ref[['GRCh38']]$pos
    tmp$BP_GRCh36<-ref[['GRCh36']]$pos
    rm(ref)

    # Insert IUPAC codes into ref
    tmp$IUPAC[tmp$A1 == 'A' & tmp$A2 =='T' | tmp$A1 == 'T' & tmp$A2 =='A']<-'W'
    tmp$IUPAC[tmp$A1 == 'C' & tmp$A2 =='G' | tmp$A1 == 'G' & tmp$A2 =='C']<-'S'
    tmp$IUPAC[tmp$A1 == 'A' & tmp$A2 =='G' | tmp$A1 == 'G' & tmp$A2 =='A']<-'R'
    tmp$IUPAC[tmp$A1 == 'C' & tmp$A2 =='T' | tmp$A1 == 'T' & tmp$A2 =='C']<-'Y'
    tmp$IUPAC[tmp$A1 == 'G' & tmp$A2 =='T' | tmp$A1 == 'T' & tmp$A2 =='G']<-'K'
    tmp$IUPAC[tmp$A1 == 'A' & tmp$A2 =='C' | tmp$A1 == 'C' & tmp$A2 =='A']<-'M'

    for(pop in unique(pop_dat$POP)){
        # Read in reference frequency data
        freq<-fread(paste0(output_dir,'/freq_files/',pop,'/ref.',pop,'.chr',chr,'.frq'))
        
        # The freq files have come from the reference files, so we can assume they are on the same strand
        freq_match<-merge(tmp, freq[,c('SNP','A1','A2','MAF'), with=F], by=c('SNP','A1','A2'))
        freq_swap<-merge(tmp, freq[,c('SNP','A1','A2','MAF'), with=F], by.x=c('SNP','A1','A2'), by.y=c('SNP','A2','A1'))
        freq_swap$MAF<-1-freq_swap$MAF
        tmp_freq<-rbind(freq_match, freq_swap)
        tmp_freq<-tmp_freq[match(tmp$SNP, tmp_freq$SNP),]
        
        tmp[[paste0('REF.FRQ.',pop)]]<-tmp_freq$MAF
    }

    tmp<-tmp[,c("CHR","SNP","BP_GRCh36","BP_GRCh37","BP_GRCh38","A1","A2","IUPAC",paste0('REF.FRQ.',unique(pop_dat$POP))), with=F]
    saveRDS(tmp, file = paste0(output_dir,'/ref.chr',chr,'.rds'))
}

# Compress folder and upload to zenodo
setwd('resources/data')
system(paste0('tar -czvf genopredpipe_1kg.tar.gz ref'))
