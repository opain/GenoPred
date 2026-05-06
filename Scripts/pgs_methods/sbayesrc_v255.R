#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NULL, type='character',
      help="Path to per chromosome reference PLINK files [required]"),
  make_option("--ref_pcs", action="store", default=NULL, type='character',
      help="Reference PCs for continuous ancestry correction [optional]"),
  make_option("--pop_data", action="store", default=NULL, type='character',
      help="File containing the population code and location of the keep file [required]"),
  make_option("--plink2", action="store", default='plink2', type='character',
      help="Path PLINK v2 software binary [required]"),
  make_option("--output", action="store", default=NULL, type='character',
      help="Path for output files [required]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
      help="Number of cores for parallel computing [optional]"),
  make_option("--sumstats", action="store", default=NULL, type='character',
      help="GWAS summary statistics in LDSC format [required]"),
  make_option("--test", action="store", default=NA, type='character',
      help="Specify number of SNPs to include [optional]"),
  make_option("--seed", action="store", default=1, type='numeric',
      help="Seed number [optional]"),
  make_option("--gctb", action="store", default=NULL, type='character',
      help="GCTB v2.5.5 binary [required]"),
  make_option("--sbayesrc_ldref", action="store", default=NULL, type='character',
      help="Path to SBayesRC LD reference data [required]"),
  make_option("--sbayesrc_annot", action="store", default=NULL, type='character',
      help="Path to SBayesRC annotations [required]")
)

opt = parse_args(OptionParser(option_list = option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')

# Check required inputs
if(is.null(opt$ref_plink_chr)){
  stop('--ref_plink_chr must be specified.\n')
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
if(is.null(opt$gctb)){
  stop('--gctb must be specified.\n')
}
if(is.null(opt$sbayesrc_ldref)){
  stop('--sbayesrc_ref must be specified.\n')
}
if(is.null(opt$sbayesrc_annot)){
  stop('--sbayesrc_annot must be specified.\n')
}

# Create output directory
opt$output_dir <- paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

# Create temp directory
tmp_dir <- tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'sbayesrc_v255.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

set.seed(opt$seed)

#####
# Read in sumstats
#####

log_add(log_file = log_file, message = 'Reading in GWAS.')

# Read in, check and format GWAS summary statistics
gwas <- read_sumstats(sumstats = opt$sumstats, chr = CHROMS, log_file = log_file, req_cols = c('SNP','A1','A2','FREQ','BETA','SE','P','N'))

###
# Change to COJO format
###

gwas <- gwas[, c('SNP','A1','A2','FREQ','BETA','SE','P','N'), with=F]
names(gwas) <- c('SNP','A1','A2','freq','b','se','p','N')

# Write out cojo format sumstats
fwrite(gwas, paste0(tmp_dir,'/GWAS_sumstats_COJO.txt'), sep=' ', na = "NA", quote=F)

rm(gwas)
gc()

# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

###
# Perform sumstat QC using SBayesRC package
###

tidy_log <- system(paste0(opt$gctb, ' ',
  '--ldm-eigen ', opt$sbayesrc_ldref,' ',
  '--gwas-summary ', tmp_dir, '/GWAS_sumstats_COJO.txt ',
  '--impute-summary ',
  '--out ',tmp_dir, '/tidy.ma ',
  '--thread ', opt$n_cores),
  intern = TRUE)

writeLines(tidy_log, paste0(opt$output,'.tidy.log'))

tidy_log_qc <- tidy_log[grepl('matched SNPs in the GWAS summary data', tidy_log)]
tidy_log_qc <- as.numeric(gsub(' .*', '', tidy_log_qc))
log_add(
  log_file = log_file,
  message = paste0(tidy_log_qc, ' variants remain after SBayesRC tidy step.')
)

# Check number of imputed variants
nsnp_imp <- system(paste0('wc -l ', tmp_dir, '/tidy.ma.imputed.ma'), intern = T)
nsnp_imp <- as.numeric(gsub(' .*', '', nsnp_imp)) - 1

log_add(
  log_file = log_file,
  message = paste0(nsnp_imp, ' variants remain after SBayesRC imputation step.')
)

###
# Run SBayesRC
###

log_add(log_file = log_file, message = 'Running SBayesRC main analysis...')

sbayerc_log <- system(paste0(
  opt$gctb, ' ', 
  '--ldm-eigen ', opt$sbayesrc_ldref,' ',
  '--gwas-summary ', tmp_dir, '/tidy.ma.imputed.ma ', 
  '--sbayes RC ',
  '--annot ', opt$sbayesrc_annot,' ',
  '--out ', tmp_dir, '/sbrc ',
#  '--chain-length 100 --burn-in 50 ', # TEST
  '--thread ', opt$n_cores), intern = T)

writeLines(sbayerc_log, paste0(opt$output,'.sbayesrc.log'))

system(paste0('cp ', tmp_dir, '/sbrc* ', opt$output_dir))

###
# Process score file
###

log_add(log_file = log_file, message = 'Processing score file...')

score <- fread(paste0(tmp_dir, '/sbrc.snpRes'))

log_add(log_file = log_file, message = paste0('Score file contains ', nrow(score), ' variants.'))
log_add(log_file = log_file, message = paste0(sum(score$A1Effect != 0, na.rm = T), ' variants have non-zero effect.'))

# Save in plink score format
score <- score[,c('Name', 'A1', 'A2', 'A1Effect'), with=F]
names(score) <- c('SNP', 'A1', 'A2', 'SCORE_SBayesRC')

# Save unmapped version of the score file
fwrite(score, paste0(opt$output,'.unmapped.score.gz'), col.names=T, sep=' ', quote=F)

# Flip effects to match reference alleles
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]
score_new <- map_score(ref = ref, score = score)

log_add(log_file = log_file, message = paste0(sum(score$SNP %in% ref$SNP), ' variants present in reference.'))
log_add(log_file = log_file, message = paste0(sum(score$SNP[score$SCORE_SBayesRC != 0] %in% ref$SNP), ' variants with non-zero effect present in reference.'))

fwrite(score_new, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

if(file.exists(paste0(opt$output,'.score.gz'))){
  system(paste0('rm ',opt$output,'.score.gz'))
}
system(paste0('gzip ',opt$output,'.score'))

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

