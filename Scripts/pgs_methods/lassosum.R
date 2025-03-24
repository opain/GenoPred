#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
	make_option("--ref_plink_chr", action="store", default=NULL, type='character',
			help="Path to genome-wide reference PLINK files [required]"),
	make_option("--ref_keep", action="store", default=NULL, type='character',
			help="Keep file to subset individuals in reference for clumping [optional]"),
	make_option("--ref_pcs", action="store", default=NULL, type='character',
	    help="Reference PCs for continuous ancestry correction [optional]"),
	make_option("--gwas_pop", action="store", default=NULL, type='character',
			help="Population of GWAS sample [required]"),
  make_option("--pop_data", action="store", default=NULL, type='character',
      help="File containing the population code and location of the keep file [required]"),
	make_option("--plink2", action="store", default='plink2', type='character',
	    help="Path PLINK v2 software binary [optional]"),
	make_option("--output", action="store", default=NULL, type='character',
			help="Path for output files [required]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
	    help="Number of cores to use [optional]"),
	make_option("--test", action="store", default=NA, type='character',
	    help="Specify number of SNPs to include [optional]"),
	make_option("--sumstats", action="store", default=NULL, type='character',
			help="GWAS summary statistics in LDSC format [optional]")
)

opt = parse_args(OptionParser(option_list = option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')
library(lassosum)
library(parallel)
cl <- makeCluster(opt$n_cores)

# Store original working directory
orig_wd<-getwd()

# Check required inputs
if(is.null(opt$ref_plink_chr)){
  stop('--ref_plink_chr must be specified.\n')
}
if(is.null(opt$gwas_pop)){
  stop('--gwas_pop must be specified.\n')
}
if(is.null(opt$sumstats)){
  stop('--sumstats must be specified.\n')
}
if(is.null(opt$pop_data)){
  stop('--pop_data must be specified.\n')
}
if(is.null(opt$output)){
  stop('--output must be specified.\n')
}

# Create output directory
opt$output_dir <- paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'lassosum.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

#####
# Read in sumstats
#####

log_add(log_file = log_file, message = 'Reading in GWAS.')

# Read in, check and format GWAS summary statistics
gwas <- read_sumstats(sumstats = opt$sumstats, chr = CHROMS, log_file = log_file, req_cols = c('CHR','SNP','BP','A1','A2','BETA','P','N'))

# Store average sample size
gwas_N <- mean(gwas$N)

###
# Merge the per chromosome reference genetic data and subset opt$ref_keep
###

log_add(log_file = log_file, message = 'Merging per chromosome reference data.')

# Save in plink1 format for lassosum
plink_merge(pfile = opt$ref_plink_chr, chr = CHROMS, plink2 = opt$plink2, keep = opt$ref_keep, extract = gwas$SNP, make_bed =T, out = paste0(tmp_dir, '/ref_merge'))

# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

#####
# Calculate correlation between SNP and phenotype
#####

cor <- p2cor(p = gwas$P, n = gwas_N, sign = gwas$BETA)

#####
# Perform lassosum to shrink effects using a range of parameters
#####

log_add(log_file = log_file, message = 'Running lassosum pipeline.')

# Change directory to lassosum package for access to LDblock data
setwd(system.file("data", package="lassosum"))

# Identify LD block data to be used
if(opt$gwas_pop %in% c('EUR','AFR')){
  ld_block_dat <- paste0(opt$gwas_pop,'.hg19')
}
if(opt$gwas_pop %in% 'EAS'){
  ld_block_dat <- 'ASN.hg19'
}
if(opt$gwas_pop %in% c('AMR','SAS')){
  ld_block_dat <- 'EUR.hg19'
  log_add(log_file = log_file, message = 'Using LD block data for EUR.')
}

# Run pipeline
out <- lassosum.pipeline(
  cor = cor,
  chr = gwas$CHR,
  pos = gwas$BP,
  A1 = gwas$A1,
  A2 = gwas$A2,
  ref.bfile = paste0(tmp_dir, '/ref_merge'),
  LDblocks = ld_block_dat,
  cluster = cl)

# Change working directory back to the original
setwd(orig_wd)

# Write out a score file
score_file <- data.table(SNP = gwas$SNP[out$sumstats$order], out$sumstats[c('A1', 'A2')])

for(i in 1:length(out$s)){
  for(k in 1:length(out$lambda)){
    score_file_tmp<-data.table(out$beta[[i]][,k])
    names(score_file_tmp)<-paste0('SCORE_s', out$s[i], '_lambda', out$lambda[k])
    score_file<-cbind(score_file, score_file_tmp)
  }
}

# Flip effects to match reference alleles
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]
score_new <- map_score(ref = ref, score = score_file)

# Reduce number of significant figures to save space
score_new[, (4:ncol(score_new)) := lapply(.SD, signif, digits = 7), .SDcols = 4:ncol(score_new)]

fwrite(score_new, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

if(file.exists(paste0(opt$output,'.score.gz'))){
  system(paste0('rm ',opt$output,'.score.gz'))
}

system(paste0('gzip ',opt$output,'.score'))

#####
# Perform pseudovalidation
#####

log_add(log_file = log_file, message = 'Performing pseudovalidation.')

# Change directory to lassosum package for access to LDblock data
setwd(system.file("data", package="lassosum"))
# Run pseudovalidation
v <- pseudovalidate(out, plot = F)
# Change working directory back to the original
setwd(orig_wd)

# Subset the validated lassosum model
out2 <- subset(out, s=v$best.s, lambda=v$best.lambda)

log_add(log_file = log_file, message = c(
  'Pseudovalidated parameters:',
  paste0('s = ', out2$s),
  paste0('lambda = ', out2$lambda),
  paste0('value = ', v$validation.table$value[v$validation.table$lambda == v$best.lambda & v$validation.table$s == v$best.s])
  ))

# Record end time of test
if(!is.na(opt$test)){
  test_finish(log_file = log_file, test_start.time = test_start.time)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()
