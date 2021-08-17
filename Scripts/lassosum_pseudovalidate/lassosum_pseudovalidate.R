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
	make_option("--plink", action="store", default='plink', type='character',
	    help="Path PLINK software binary [required]"),
	make_option("--output", action="store", default='./Output', type='character',
			help="Path for output files [required]"),
	make_option("--test", action="store", default=NA, type='character',
	    help="Specify number of SNPs to include [optional]"),
	make_option("--n_cores", action="store", default=1, type='numeric',
	    help="Number of cores to use [optional]"),
	make_option("--prune_mhc", action="store", default=F, type='logical',
	    help="Set to T if MHC region should be pruned [optional]"),
	make_option("--sumstats", action="store", default=NA, type='character',
			help="GWAS summary statistics [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(lassosum)
library(parallel)
cl <- makeCluster(opt$n_cores)
orig_wd<-getwd()

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)
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
# lassosum_pseudovalidate.R V1.0
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

# Remove MHC
if(opt$prune_mhc == T){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Extracted top variant in MHC region.\n')
  sink()
  
  GWAS_hla<-GWAS[GWAS$CHR == 6 & GWAS$BP > 28e6 & GWAS$BP < 34e6,]
  GWAS_hla_excl<-GWAS_hla$SNP[GWAS_hla$P != min(GWAS_hla$P)]
  GWAS<-GWAS[!(GWAS$SNP %in% GWAS_hla_excl),]
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
                       LDblocks = 'EUR.hg19', cluster=cl)

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
v <- pseudovalidate(out, cluster=cl)
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

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('value =', v$validation.table$value[v$validation.table$lambda == v$best.lambda & v$validation.table$s == v$best.s],'\n')
sink()

if(!is.na(opt$ref_keep)){
  system(paste0('rm ',opt$output_dir,'lassosum_ref_gw'))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
