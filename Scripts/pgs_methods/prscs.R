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
make_option("--PRScs_path", action="store", default=NA, type='character',
		help="Path to PRScs executable [required]"),
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Number of cores for parallel computing [optional]"),
make_option("--PRScs_ref_path", action="store", default=T, type='character',
		help="Path to PRScs reference [required]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--phi_param", action="store", default='auto', type='character',
    help="Path to PRScs reference [optional]"),
make_option("--seed", action="store", default=NA, type='numeric',
    help="Seed number for PRScs [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(foreach)
library(doMC)
registerDoMC(opt$n_cores)
source('../Scripts/functions/misc.R')

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

phi_param<-unlist(strsplit(opt$phi_param,','))

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

if(any(grepl('auto',phi_param))){
	phi_param<-c(sprintf("%1.00e", as.numeric(phi_param[!grepl('auto',phi_param)])),'auto')
} else {
	phi_param<-sprintf("%1.00e", as.numeric(phi_param))
}

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# prscs.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#####
# Read in sumstats and insert p-values
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in GWAS and harmonising with reference.\n')
sink()

GWAS<-fread(cmd=paste0('zcat ',opt$sumstats), nThread=1)
GWAS<-GWAS[complete.cases(GWAS),]

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

GWAS_N<-mean(GWAS$N)

if(('BETA' %in% names(GWAS))){
  GWAS<-GWAS[,c('SNP','A1','A2','BETA','P')]
} else {
  GWAS<-GWAS[,c('SNP','A1','A2','OR','P')]
}

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

#####
# Process sumstats using PRSsc
#####

# Make a data.frame listing chromosome and phi combinations
jobs<-NULL
for(i in CHROMS){
  jobs<-rbind(jobs,data.frame(CHR=i,
                       phi=phi_param))
}

# Run using PRScs auto, and specifying a range of global shrinkage parameters (1e-6, 1e-4, 1e-2, 1)
foreach(i=1:dim(jobs)[1], .combine=c, .options.multicore=list(preschedule=FALSE)) %dopar% {
  if(is.na(opt$seed)){
    if(jobs$phi[i] == 'auto'){
      system(paste0(opt$PRScs_path,' --ref_dir=',opt$PRScs_ref_path,' --bim_prefix=',opt$ref_plink_chr,jobs$CHR[i],' --sst_file=',opt$output_dir,'GWAS_sumstats_temp.txt --n_gwas=',round(GWAS_N,0),' --out_dir=',opt$output,' --chrom=',jobs$CHR[i]))
    } else {
      system(paste0(opt$PRScs_path,' --ref_dir=',opt$PRScs_ref_path,' --bim_prefix=',opt$ref_plink_chr,jobs$CHR[i],' --phi=',jobs$phi[i],' --sst_file=',opt$output_dir,'GWAS_sumstats_temp.txt --n_gwas=',round(GWAS_N,0),' --out_dir=',opt$output,' --chrom=',jobs$CHR[i]))
    }
  } else {
    if(jobs$phi[i] == 'auto'){
      system(paste0(opt$PRScs_path,' --ref_dir=',opt$PRScs_ref_path,' --bim_prefix=',opt$ref_plink_chr,jobs$CHR[i],' --sst_file=',opt$output_dir,'GWAS_sumstats_temp.txt --n_gwas=',round(GWAS_N,0),' --out_dir=',opt$output,' --chrom=',jobs$CHR[i],' --seed=',opt$seed))
    } else {
      system(paste0(opt$PRScs_path,' --ref_dir=',opt$PRScs_ref_path,' --bim_prefix=',opt$ref_plink_chr,jobs$CHR[i],' --phi=',jobs$phi[i],' --sst_file=',opt$output_dir,'GWAS_sumstats_temp.txt --n_gwas=',round(GWAS_N,0),' --out_dir=',opt$output,' --chrom=',jobs$CHR[i],' --seed=',opt$seed))
    }
  }
}

system(paste0('rm ',opt$output_dir,'GWAS_sumstats_temp.txt'))

if(!is.na(opt$test)){
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
}

####
# Combine score files 
####
score_all<-NULL
for(phi_i in phi_param){
  score_phi<-NULL
  for(i in CHROMS){
    score_phi_i<-fread(paste0(opt$output,'_pst_eff_a1_b0.5_phi',phi_i,'_chr',i,'.txt'))
    score_phi<-rbind(score_phi, score_phi_i)
  }
  if(phi_i == phi_param[1]){
    score_phi<-score_phi[,c('V2','V4','V6'), with=F]
    names(score_phi)<-c('SNP','A1',paste0('SCORE_phi_',phi_i))
  } else {
    score_phi<-score_phi[,'V6', with=F]
    names(score_phi)<-paste0('SCORE_phi_',phi_i)
  }
  score_all<-cbind(score_all, score_phi)
}

fwrite(score_all, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

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

system(paste0('rm ',opt$output,'_pst_eff_a1_b0.5_*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
