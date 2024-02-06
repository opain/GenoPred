#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NA, type='character',
              help="Path to per chromosome reference PLINK files [required]"),
  make_option("--ref_pop_scale", action="store", default=NA, type='character',
              help="File containing the population code and location of the keep file [required]"),
  make_option("--plink2", action="store", default='plink', type='character',
              help="Path PLINK v2 software binary [required]"),
  make_option("--output", action="store", default='NA', type='character',
              help="Path for output files [required]"),
  make_option("--memory", action="store", default=5000, type='numeric',
              help="Memory limit [optional]"),
  make_option("--sumstats", action="store", default=NA, type='character',
              help="GWAS summary statistics [optional]"),
  make_option("--ldak", action="store", default=NA, type='character',
              help="Path to ldak executable [required]"),
  make_option("--quick_prs_ref", action="store", default=NA, type='character',
              help="Path to folder containing ldak quick prs reference [required]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
              help="Number of cores for parallel computing [optional]"),
  make_option("--prs_model", action="store", default='bayesr', type='character',
              help="Model used for deriving SNP-weights [optional]"),
  make_option("--genomic_control", action="store", default=F, type='logical',
              help="Logical indicating whether genomic control was applied to GWAS [optional]"),
  make_option("--test", action="store", default=NA, type='character',
              help="Specify number of SNPs to include [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

CHROMS<-1:22

if(!is.na(opt$test)){
  if(grepl('chr', opt$test)){
    single_chr_test<-T
    CHROMS<-as.numeric(gsub('chr','',opt$test))
  } else {
    single_chr_test<-F
    opt$test<-as.numeric(opt$test)
  }
}

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
  '#################################################################
# quickprs.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#####
# Format the sumstats
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in GWAS and harmonising with reference.\n')
sink()

GWAS<-fread(cmd=paste0('zcat ',opt$sumstats), nThread=opt$n_cores)
GWAS<-GWAS[complete.cases(GWAS),]

# Update BP to match the reference
ref_bim<-NULL
for(i in 1:22){
  ref_bim<-rbind(ref_bim, fread(paste0(opt$ref_plink_chr,i,'.bim'), header=F))
}

GWAS$BP<-NULL
GWAS<-merge(GWAS, ref_bim[,c('V2','V4'), with=F], by.x='SNP',by.y='V2')
names(GWAS)[names(GWAS) == 'V4']<-'BP'
GWAS<-GWAS[order(GWAS$CHR, GWAS$BP),]

# Extract subset if testing
if(!is.na(opt$test)){
  if(single_chr_test == F){
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Testing mode enabled. Extracted ',opt$test,' variants per chromsome.\n', sep='')
    sink()
    
    GWAS_test<-NULL
    for(i in 1:22){
      GWAS_tmp<-GWAS[GWAS$CHR == i,]
      GWAS_tmp<-GWAS_tmp[order(GWAS_tmp$BP),]
      GWAS_tmp<-GWAS_tmp[1:opt$test,]
      GWAS_test<-rbind(GWAS_test,GWAS_tmp)
    }
    
    GWAS<-GWAS_test
    GWAS<-GWAS[complete.cases(GWAS),]
    rm(GWAS_test)
    print(table(GWAS$CHR))
    
  } else {
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Testing mode enabled. Extracted chromosome ',opt$test,' variants per chromsome.\n', sep='')
    sink()
    
    GWAS<-GWAS[GWAS$CHR == CHROMS,]
    print(table(GWAS$CHR))
  }
}

# Update RSID to CHR:BP to match reference data
# During sumstat cleaning, these columns are replaced with those found in build hg19 which is appropriate
GWAS$SNP<-paste0(GWAS$CHR,':',GWAS$BP)

# Check overlap with LDAK HapMap3 SNP-list
ldak_hm3_file<-list.files(opt$quick_prs_ref)
ldak_hm3_file<-ldak_hm3_file[grepl('.cors.bim',ldak_hm3_file)]
ldak_hm3<-fread(paste0(opt$quick_prs_ref,'/',ldak_hm3_file))

if(is.na(opt$test)){
  ref_overlap<-sum(GWAS$SNP %in% ldak_hm3$V2)/nrow(ldak_hm3)
  GWAS<-GWAS[GWAS$SNP %in% ldak_hm3$V2,]
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('GWAS-reference overlap is ',round(ref_overlap*100,2),'%.\n', sep='')
  sink()
} else {
  GWAS<-GWAS[GWAS$SNP %in% ldak_hm3$V2,]
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('GWAS-reference overlap check is skipped as this is a test run.\n', sep='')
  sink()
}

# Calculate Z
if(('BETA' %in% names(GWAS))){
  GWAS$Z<-abs(GWAS$BETA)/GWAS$SE
  GWAS$Z[GWAS$BETA < 0]<- -GWAS$Z[GWAS$BETA < 0]
} else {
  GWAS$Z<-abs(log(GWAS$OR))/GWAS$SE
  GWAS$Z[GWAS$OR < 1]<- -GWAS$Z[GWAS$OR < 1]
}

GWAS$Predictor<-paste0(GWAS$CHR, ':', GWAS$BP)

GWAS<-GWAS[,c('Predictor','A1','A2','N','Z')]
names(GWAS)<-c('Predictor','A1','A2','n','Z')

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

fwrite(GWAS, paste0(opt$output_dir,'GWAS_sumstats_temp.txt'), sep=' ')

rm(GWAS)
gc()

if(!is.na(opt$test)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  test_start.time <- Sys.time()
  cat('Test started at',as.character(test_start.time),'\n')
  sink()
}

############
# Estimate Per-Predictor Heritabilities
############

# Calculate Per-Predictor Heritabilities.
ref_files<-list.files(opt$quick_prs_ref)

tagging_file<-ref_files[grepl('quickprs.tagging',ref_files)]
matrix_file<-ref_files[grepl('quickprs.matrix',ref_files)]

if(opt$genomic_control == F){
  system(paste0(opt$ldak,' --sum-hers ',opt$output_dir,'/bld.ldak --tagfile ',opt$quick_prs_ref,'/',tagging_file,' --summary ',opt$output_dir,'GWAS_sumstats_temp.txt --matrix ',opt$quick_prs_ref,'/',matrix_file,' --max-threads ',opt$n_cores,' --check-sums NO'))
} else{
  system(paste0(opt$ldak,' --sum-hers ',opt$output_dir,'/bld.ldak --genomic-control YES --tagfile ',opt$quick_prs_ref,'/',tagging_file,' --summary ',opt$output_dir,'GWAS_sumstats_temp.txt --matrix ',opt$quick_prs_ref,'/',matrix_file,' --max-threads ',opt$n_cores,' --check-sums NO'))
}

ldak_res_her<-fread(paste0(opt$output_dir,'/bld.ldak.hers'))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('SNP-based heritability estimated to be ',ldak_res_her$Heritability[nrow(ldak_res_her)]," (SD=", ldak_res_her$SD[nrow(ldak_res_her)],").\n",sep='')
sink()

######
# Estimate effect sizes for training and full prediction models.
######

cor_file_prefix<-gsub('.cors.bin','',ref_files[grepl('.cors.bin',ref_files)])

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Running MegaPRS: ',opt$prs_model,' model.\n', sep='')
sink()

system(paste0(opt$ldak,' --mega-prs ',opt$output_dir,'/mega_full --model ',opt$prs_model,' --cors ',opt$quick_prs_ref,'/',cor_file_prefix,' --ind-hers ',opt$output_dir,'/bld.ldak.ind.hers --summary ',opt$output_dir,'GWAS_sumstats_temp.txt --high-LD ',opt$quick_prs_ref,'/highld.snps --cv-proportion 0.1 --window-cm 1 --max-threads ',opt$n_cores,' --extract ',opt$output_dir,'GWAS_sumstats_temp.txt'))

# Save the parameters file
system(paste0('cp ',opt$output_dir,'/mega_full.parameters ',opt$output,'.model_param.txt'))

# Save the pseudosummary results
system(paste0('cp ',opt$output_dir,'/mega_full.cors ',opt$output,'.pseudoval.txt'))

# Identify the best fitting model
ldak_res_cors<-fread(paste0(opt$output_dir,'/mega_full.cors'), nThread=opt$n_cores)
best_score<-ldak_res_cors[ldak_res_cors$Correlation == max(ldak_res_cors$Correlation),]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Model/Models ',paste0(best_score$Model, collapse=', '),' was/were identified as the best with correlation of ',best_score$Correlation[1],'.\n', sep='')
sink()

if(!is.na(opt$test)){
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
}

######
# Format final score file
######

# Read in the scores
score<-fread(paste0(opt$output_dir,'/mega_full.effects'), nThread=opt$n_cores)

# Change IDs to RSIDs
bim<-NULL
for(i in 1:22){
  bim<-rbind(bim, fread(paste0(opt$ref_plink_chr,i,'.bim'), nThread=opt$n_cores, header=F))
}
bim$Predictor<-paste0(bim$V1,':',bim$V4)
score<-merge(score, bim[,c('Predictor','V2'),with=F], by='Predictor')
score<-score[,c('V2','A1',names(score)[grepl('Model', names(score))]), with=F]
names(score)[1]<-'SNP'
names(score)[grepl('Model', names(score))]<-paste0('SCORE_ldak_',names(score)[grepl('Model', names(score))])

fwrite(score, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

if(file.exists(paste0(opt$output,'.score.gz'))){
  system(paste0('rm ',opt$output,'.score.gz'))
}

system(paste0('gzip ',opt$output,'.score'))

####
# Calculate mean and sd of polygenic scores
####

# Calculate polygenic scores for reference individuals
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

scores<-calc_score(
  bfile=opt$ref_plink_chr, 
  score=paste0(opt$output,'.score.gz')
)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

# Calculate the mean and sd of scores for each population specified in pop_scale
pop_keep_files<-read.table(opt$ref_pop_scale, header=F, stringsAsFactors=F)

for(k in 1:dim(pop_keep_files)[1]){
	pop<-pop_keep_files$V1[k]
	keep<-fread(pop_keep_files$V2[k], header=F)
  names(keep)<-c('FID','IID')
  ref_scale<-score_mean_sd(scores=scores, keep=keep)
	fwrite(ref_scale, paste0(opt$output,'.',pop,'.scale'), sep=' ')
}

###
# Clean up temporary files
###

system(paste0('rm ',opt$output_dir,'bld*'))
system(paste0('rm ',opt$output_dir,'GWAS_sums*'))
system(paste0('rm ',opt$output_dir,'mega*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()

