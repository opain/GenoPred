#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
	make_option("--ref_plink_gw", action="store", default=NA, type='character',
			help="Path to genome-wide reference PLINK files [required]"),
	make_option("--ref_keep", action="store", default=NA, type='character',
			help="Keep file to subset individuals in reference for clumping [required]"),
	make_option("--ref_freq_chr", action="store", default=NA, type='character',
	    help="Path to per chromosome PLINK frequency (.frq) files [optional]"),
	make_option("--ref_maf", action="store", default=NA, type='numeric',
	    help="Minor allele frequency threshold to be applied based on ref_freq_chr [optional]"),
	make_option("--ref_pop_scale", action="store", default=NA, type='character',
			help="File containing the population code and location of the keep file [required]"),
	make_option("--plink2", action="store", default='plink2', type='character',
	    help="Path PLINK v2 software binary [required]"),
	make_option("--output", action="store", default='./Output', type='character',
			help="Path for output files [required]"),
	make_option("--memory", action="store", default=5000, type='numeric',
	    help="Memory limit [optional]"),
	make_option("--test", action="store", default=NA, type='character',
	    help="Specify number of SNPs to include [optional]"),
	make_option("--sumstats", action="store", default=NA, type='character',
			help="GWAS summary statistics in LDSC format [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(lassosum)
orig_wd<-getwd()

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
# polygenic_score_file_creator_lassosum.R V1.0
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
GWAS<-GWAS[GWAS$P > 0,]
GWAS_N<-mean(GWAS$N)

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

if(!('BETA' %in% names(GWAS))){
  GWAS$BETA<-log(GWAS$OR)
}

GWAS<-GWAS[,c('CHR','SNP','BP','A1','A2','BETA','P')]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

if(!is.na(opt$ref_keep)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('ref_keep used to subset reference genotype data.\n')
  sink()
  
  #####
  # Create subset of ref files
  #####
  
  system(paste0(opt$plink,' --bfile ',opt$ref_plink_gw,' --keep ',opt$ref_keep,' --make-bed --out ',opt$output_dir,'lassosum_ref_gw'))
  
  opt$ref_plink_subset<-paste0(opt$output_dir,'lassosum_ref_gw')
} else {
  opt$ref_plink_subset<-opt$ref_plink_gw
}

if(!is.na(opt$test)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  test_start.time <- Sys.time()
  cat('Test started at',as.character(test_start.time),'\n')
  sink()
}

#####
# Calculate correlation between SNP and phenotype 
#####

cor <- p2cor(p = GWAS$P, n = GWAS_N, sign=GWAS$BETA)

#####
# Perform lassosum to shrink effects using a range of parameters
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Running lassosum pipeline...')
sink()

setwd(system.file("data", package="lassosum"))

out<-lassosum.pipeline(cor=cor, chr=GWAS$CHR, pos=GWAS$BP, 
                       A1=GWAS$A1, A2=GWAS$A2,
                       ref.bfile=paste0(orig_wd,'/',opt$ref_plink_subset), 
                       LDblocks = 'EUR.hg19')

setwd(orig_wd)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

rm(cor)
gc()

#####
# Perform pseudovalidation to idenitfy the best p-value threshold
#####
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Idenitfying best parameters via pseudovalidation...')
sink()

bitmap(paste0(opt$output,'.pseudovalidate.png'), unit='px', res=300, height=2000, width=2000)
setwd(system.file("data", package="lassosum"))
v <- pseudovalidate(out)
setwd(orig_wd)
dev.off()

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

# Subset the validated lassosum model
out2 <- subset(out, s=v$best.s, lambda=v$best.lambda)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Pseudovalidated parameters:
s = ',out2$s,'
lambda = ',out2$lambda,'\n',sep='')
sink()

rm(v)
gc()

# Write out a score file
score_file<-data.table(SNP=GWAS$SNP[out$sumstats$order],
                       out$sumstats[c('A1')])
score_file<-score_file[,c('SNP','A1'),with=F]

for(i in 1:length(out$s)){
  for(k in 1:length(out$lambda)){
    score_file_tmp<-data.table(out$beta[[i]][,k])
    names(score_file_tmp)<-paste0('SCORE_s',out$s[i],'_lambda',out$lambda[k])
    score_file<-cbind(score_file,score_file_tmp)
  }
}

fwrite(score_file, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

if(file.exists(paste0(opt$output,'.score.gz'))){
  system(paste0('rm ',opt$output,'.score.gz'))
}

system(paste0('gzip ',opt$output,'.score'))

if(!is.na(opt$test)){
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  system(paste0('rm ',opt$output,'*.score'))
  system(paste0('rm ',opt$output,'.pseudovalidate.png'))
  q()
}

#####
# Calculate the mean and sd of scores for each population specified in pop_scale
#####

# Calculate polygenic scores for reference individuals
# Do this using PLINK not lassosum as we want to use MAF to impute missing genotypes in target samples, which lassosum does not do.
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

if(ncol(score_file)-2 == 1){
  for(i in CHROMS){
    system(paste0(opt$plink2, ' --bfile ',opt$ref_plink_gw,' --score ',opt$output,'.score.gz header-read cols=fid,nallele,denom,dosagesum,scoresums --out ',opt$output,'.profiles.GW --memory ',floor(opt$memory*0.7)))
  }
} else {
  for(i in CHROMS){
    system(paste0(opt$plink2, ' --bfile ',opt$ref_plink_gw,' --score ',opt$output,'.score.gz header-read cols=fid,nallele,denom,dosagesum,scoresums --score-col-nums 3-',2+ncol(score_file)-2,' --out ',opt$output,'.profiles.GW --memory ',floor(opt$memory*0.7)))
  }
}

sscore<-fread(paste0(opt$output,'.profiles.GW.sscore'))
fam<-sscore[,1:2,with=F]
scores<-sscore[,grepl('SCORE_.*_SUM', names(sscore)),with=F]
scores<-cbind(fam,scores)

names(scores)<-c('FID','IID',names(score_file)[-1:-2])

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

#####
# Calculate the mean and sd of scores for each population specified in pop_scale
#####
pop_keep_files<-read.table(opt$ref_pop_scale, header=F, stringsAsFactors=F)

for(k in 1:dim(pop_keep_files)[1]){
	pop<-pop_keep_files$V1[k]
	keep<-fread(pop_keep_files$V2[k], header=F)
	scores_keep<-scores[(scores$FID %in% keep$V1),]

	ref_scale<-data.frame(	Param=names(scores_keep[,-1:-2]),
													Mean=round(sapply(scores_keep[,-1:-2], function(x) mean(x)),3),
													SD=round(sapply(scores_keep[,-1:-2], function(x) sd(x)),3))

	fwrite(ref_scale, paste0(opt$output,'.',pop,'.scale'), sep=' ')
}

if(!is.na(opt$ref_keep)){
  system(paste0('rm ',opt$output_dir,'lassosum_ref_gw*'))
}
system(paste0('rm ',opt$output,'.profiles*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
