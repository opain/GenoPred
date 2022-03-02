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
make_option("--sumstats1", action="store", default=NA, type='character',
    help="GWAS summary statistics 1 [required]"),
make_option("--sumstats2", action="store", default=NA, type='character',
    help="GWAS summary statistics 2 [required]"),
make_option("--pop1", action="store", default=NA, type='character',
    help="Super population or  sumstats1 [required]"),
make_option("--pop2", action="store", default=NA, type='character',
    help="Super population for sumstats2 [required]"),
make_option("--PRS_CSx_path", action="store", default=NA, type='character',
    help="Path to PRScs executable [required]"),
make_option("--PRS_CSx_ref_path", action="store", default=NA, type='character',
    help="Path to PRScs ld reference data [required]"),
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Number of cores for parallel computing [optional]"),
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

opt$output_dir<-paste0(dirname(opt$output),'/')
opt$output_name<-basename(opt$output)
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
# polygenic_score_file_creator_PRScs.R V1.0
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

###
# sumstats1
###

GWAS1<-fread(cmd=paste0('zcat ',opt$sumstats1), nThread=1)
GWAS1<-GWAS1[complete.cases(GWAS1),]

# Extract subset if testing
if(!is.na(opt$test)){
  if(single_chr_test == F){
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Testing mode enabled. Extracted ',opt$test,' variants per chromsome.\n', sep='')
    sink()
    
    GWAS1_test<-NULL
    for(i in 1:22){
      GWAS1_tmp<-GWAS1[GWAS1$CHR == i,]
      GWAS1_tmp<-GWAS1_tmp[order(GWAS1_tmp$BP),]
      GWAS1_tmp<-GWAS1_tmp[1:opt$test,]
      GWAS1_test<-rbind(GWAS1_test,GWAS1_tmp)
    }
    
    GWAS1<-GWAS1_test
    GWAS1<-GWAS1[complete.cases(GWAS1),]
    rm(GWAS1_test)
    print(table(GWAS1$CHR))
    
  } else {
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Testing mode enabled. Extracted chromosome ',opt$test,' variants per chromsome.\n', sep='')
    sink()
    
    GWAS1<-GWAS1[GWAS1$CHR == CHROMS,]
    print(table(GWAS1$CHR))
  }
}

GWAS1_N<-mean(GWAS1$N)

if(('BETA' %in% names(GWAS1))){
  GWAS1<-GWAS1[,c('SNP','A1','A2','BETA','P')]
} else {
  GWAS1<-GWAS1[,c('SNP','A1','A2','OR','P')]
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('sumstats1 contains',dim(GWAS1)[1],'variants.\n')
sink()

fwrite(GWAS1, paste0(opt$output_dir,'GWAS1_sumstats_temp.txt'), sep=' ')

rm(GWAS1)
gc()

###
# sumstats2
###

GWAS2<-fread(cmd=paste0('zcat ',opt$sumstats2), nThread=1)
GWAS2<-GWAS2[complete.cases(GWAS2),]

# Extract subset if testing
if(!is.na(opt$test)){
  if(single_chr_test == F){
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Testing mode enabled. Extracted ',opt$test,' variants per chromsome.\n', sep='')
    sink()
    
    GWAS2_test<-NULL
    for(i in 1:22){
      GWAS2_tmp<-GWAS2[GWAS2$CHR == i,]
      GWAS2_tmp<-GWAS2_tmp[order(GWAS2_tmp$BP),]
      GWAS2_tmp<-GWAS2_tmp[1:opt$test,]
      GWAS2_test<-rbind(GWAS2_test,GWAS2_tmp)
    }
    
    GWAS2<-GWAS2_test
    GWAS2<-GWAS2[complete.cases(GWAS2),]
    rm(GWAS2_test)
    print(table(GWAS2$CHR))
    
  } else {
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Testing mode enabled. Extracted chromosome ',opt$test,' variants per chromsome.\n', sep='')
    sink()
    
    GWAS2<-GWAS2[GWAS2$CHR == CHROMS,]
    print(table(GWAS2$CHR))
  }
}

GWAS2_N<-mean(GWAS2$N)

if(('BETA' %in% names(GWAS2))){
  GWAS2<-GWAS2[,c('SNP','A1','A2','BETA','P')]
} else {
  GWAS2<-GWAS2[,c('SNP','A1','A2','OR','P')]
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('sumstats1 contains',dim(GWAS2)[1],'variants.\n')
sink()

fwrite(GWAS2, paste0(opt$output_dir,'GWAS2_sumstats_temp.txt'), sep=' ')

rm(GWAS2)
gc()

if(!is.na(opt$test)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  test_start.time <- Sys.time()
  cat('Test started at',as.character(test_start.time),'\n')
  sink()
}

#####
# Process sumstats using PRScsx
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
      system(paste0(opt$PRS_CSx_path,' --meta=True --ref_dir=',opt$PRS_CSx_ref_path,' --bim_prefix=',opt$ref_plink_chr,jobs$CHR[i],' --sst_file=',opt$output_dir,'GWAS1_sumstats_temp.txt,',opt$output_dir,'GWAS2_sumstats_temp.txt --pop=',opt$pop1,',',opt$pop2,' --n_gwas=',round(GWAS1_N,0),',',round(GWAS2_N,0),' --out_dir=',opt$output_dir,' --out_name=',opt$output_name,' --chrom=',jobs$CHR[i]))
    } else {
      system(paste0(opt$PRS_CSx_path,' --meta=True --ref_dir=',opt$PRS_CSx_ref_path,' --bim_prefix=',opt$ref_plink_chr,jobs$CHR[i],' --phi=',jobs$phi[i],' --sst_file=',opt$output_dir,'GWAS1_sumstats_temp.txt,',opt$output_dir,'GWAS2_sumstats_temp.txt --pop=',opt$pop1,',',opt$pop2,' --n_gwas=',round(GWAS1_N,0),',',round(GWAS2_N,0),' --out_dir=',opt$output_dir,' --out_name=',opt$output_name,' --chrom=',jobs$CHR[i]))
    }
  } else {
    if(jobs$phi[i] == 'auto'){
      system(paste0(opt$PRS_CSx_path,' --meta=True --ref_dir=',opt$PRS_CSx_ref_path,' --bim_prefix=',opt$ref_plink_chr,jobs$CHR[i],' --sst_file=',opt$output_dir,'GWAS1_sumstats_temp.txt,',opt$output_dir,'GWAS2_sumstats_temp.txt --pop=',opt$pop1,',',opt$pop2,' --n_gwas=',round(GWAS1_N,0),',',round(GWAS2_N,0),' --out_dir=',opt$output_dir,' --out_name=',opt$output_name,' --chrom=',jobs$CHR[i],' --seed=',opt$seed))
    } else {
      system(paste0(opt$PRS_CSx_path,' --meta=True --ref_dir=',opt$PRS_CSx_ref_path,' --bim_prefix=',opt$ref_plink_chr,jobs$CHR[i],' --phi=',jobs$phi[i],' --sst_file=',opt$output_dir,'GWAS1_sumstats_temp.txt,',opt$output_dir,'GWAS2_sumstats_temp.txt --pop=',opt$pop1,',',opt$pop2,' --n_gwas=',round(GWAS1_N,0),',',round(GWAS2_N,0),' --out_dir=',opt$output_dir,' --out_name=',opt$output_name,' --chrom=',jobs$CHR[i],' --seed=',opt$seed))
    }
  }
}

system(paste0('rm ',opt$output_dir,'GWAS1_sumstats_temp.txt'))
system(paste0('rm ',opt$output_dir,'GWAS2_sumstats_temp.txt'))

if(!is.na(opt$test)){
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  system(paste0('rm ',opt$output,'*.txt'))
  q()
}

####
# Combine score files 
####
for(pop_disc in c(opt$pop1, opt$pop2, 'META')){
  score_all<-NULL
  for(phi_i in phi_param){
    score_phi<-NULL
    for(i in CHROMS){
      score_phi_i<-fread(paste0(opt$output,'_',pop_disc,'_pst_eff_a1_b0.5_phi',phi_i,'_chr',i,'.txt'))
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
  
  fwrite(score_all, paste0(opt$output,'.',pop_disc,'.score'), col.names=T, sep=' ', quote=F)
  
  if(file.exists(paste0(opt$output,'.',pop_disc,'.score.gz'))){
    system(paste0('rm ',opt$output,'.',pop_disc,'.score.gz'))
  }
  
  system(paste0('gzip ',opt$output,'.',pop_disc,'.score'))
}
  
####
# Calculate mean and sd of polygenic scores
####

# Calculate polygenic scores for reference individuals
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

for(pop_disc in c(opt$pop1, opt$pop2, 'META')){
  if(length(phi_param) == 1){
    for(i in CHROMS){
      system(paste0(opt$plink2, ' --bfile ',opt$ref_plink_chr,i,' --score ',opt$output,'.',pop_disc,'.score.gz header-read --out ',opt$output,'.',pop_disc,'.profiles.chr',i,' --memory ',floor(opt$memory*0.7)))
    }
  } else {
    for(i in CHROMS){
      system(paste0(opt$plink2, ' --bfile ',opt$ref_plink_chr,i,' --score ',opt$output,'.',pop_disc,'.score.gz header-read --score-col-nums 3-',2+length(phi_param),' --out ',opt$output,'.',pop_disc,'.profiles.chr',i,' --memory ',floor(opt$memory*0.7)))
    }
  }
}

# Add up the scores across chromosomes
fam<-fread(paste0(opt$ref_plink_chr,'22.fam'))

for(pop_disc in c(opt$pop1, opt$pop2, 'META')){
  scores<-list()
  for(i in as.character(CHROMS)){
    sscore<-fread(paste0(opt$output,'.',pop_disc,'.profiles.chr',i,'.sscore'))
    scores[[i]]<-sscore[,grepl('SCORE_', names(sscore)),with=F]
    scores[[i]]<-as.matrix(scores[[i]]*sscore$NMISS_ALLELE_CT)
  }
  
  scores<-Reduce(`+`, scores)
  scores<-data.table(FID=fam$V1,
                     IID=fam$V2,
                     scores)
  
  names(scores)<-c('FID','IID',names(score_all)[-1:-2])
  
  # Calculate the mean and sd of scores for each population specified in pop_scale
  pop_keep_files<-read.table(opt$ref_pop_scale, header=F, stringsAsFactors=F)
  
  for(k in 1:dim(pop_keep_files)[1]){
  	pop<-pop_keep_files$V1[k]
  	keep<-fread(pop_keep_files$V2[k], header=F)
  	scores_keep<-scores[(scores$FID %in% keep$V1),]
  
  	ref_scale<-data.frame(	Param=names(scores_keep[,-1:-2]),
  													Mean=round(sapply(scores_keep[,-1:-2], function(x) mean(x)),3),
  													SD=round(sapply(scores_keep[,-1:-2], function(x) sd(x)),3))
  
  	fwrite(ref_scale, paste0(opt$output,'.',pop_disc,'.',pop,'.scale'), sep=' ')
  }
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Clean up temporary files
###

system(paste0('rm ',opt$output,'*.profiles.*'))
system(paste0('rm ',opt$output,'*_pst_eff_a1_b0.5_*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
