#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NULL, type='character',
      help="Path to per chromosome reference PLINK files [required]"),
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
  make_option("--gctb", action="store", default=NULL, type='character',
      help="Path to GCTB binary [required]"),
  make_option("--test", action="store", default=NA, type='character',
      help="Specify number of SNPs to include [optional]"),
  make_option("--seed", action="store", default=1, type='numeric',
      help="Seed number [optional]"),
  make_option("--sbayesrc_ldref", action="store", default=NULL, type='character',
      help="Path to SBayesRC LD reference data [required]"),
  make_option("--sbayesrc_annot", action="store", default=NULL, type='character',
      help="Path to SBayesRC annotations [required]")
)

opt = parse_args(OptionParser(option_list = option_list))

# Load dependencies
library(SBayesRC)
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
log_header(log_file = log_file, opt = opt, script = 'sbayesrc.R', start.time = start.time)

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

SBayesRC::tidy(
  mafile = paste0(tmp_dir, '/GWAS_sumstats_COJO.txt'),
  LDdir = opt$sbayesrc_ldref,
  output = paste0(tmp_dir, '/tidy.ma'),
  log2file = TRUE
)
tidy_log <- readLines(paste0(tmp_dir, '/tidy.ma.log'))
tidy_log <- tidy_log[grepl('SNPs remained after QC', tidy_log)]
tidy_log <- as.numeric(gsub(' .*', '', tidy_log))
log_add(
  log_file = log_file,
  message = paste0(tidy_log, ' variants remain after SBayesRC tidy step.')
)

system(paste0('cp ', tmp_dir, '/tidy.ma.log ', opt$output,'.tidy.log'))

###
# Impute the GWAS sumstats using SBayesRC
###

log_add(log_file = log_file, message = 'Running SBayesRC imputation step...')

SBayesRC::impute(
  mafile = paste0(tmp_dir, '/tidy.ma'),
  LDdir = opt$sbayesrc_ldref,
  output = paste0(tmp_dir, '/tidy.imp.ma'),
  log2file = T
)

###
# Run SBayesRC
###

log_add(log_file = log_file, message = 'Running SBayesRC main analysis...')

tryCatch(
  {
    # First attempt to run sbayesrc
    SBayesRC::sbayesrc(
      mafile = paste0(tmp_dir, '/tidy.imp.ma'),
      LDdir = opt$sbayesrc_ldref,
      outPrefix = paste0(tmp_dir, '/sbrc'),
      annot = opt$sbayesrc_annot,
      log2file = FALSE
    )
  },
  error = function(e) {
    # Check if the error message matches the specific issue
    if (grepl("All correlations are negative, this may indicate errors in summary data.", e$message)) {
      message("Specific error encountered. Retrying with btune = FALSE...")
      log_add(log_file = log_file, message = 'Specific error encountered. Retrying with btune = FALSE...')
      # Retry with the btune parameter set to FALSE
      SBayesRC::sbayesrc(
        mafile = paste0(tmp_dir, '/tidy.imp.ma'),
        LDdir = opt$sbayesrc_ldref,
        outPrefix = paste0(tmp_dir, '/sbrc'),
        annot = opt$sbayesrc_annot,
        bTune = FALSE,
        log2file = FALSE
      )
    } else if (grepl("An unexpected error occurred: Warning, the best parameter is the minimumn threshold, we suggest to expand the tuning grid by specify lower tuning value", e$message)) {
        message("Specific error encountered. Retrying with tuneStep=c(0.995, 0.9, 0.8, 0.7, 0.6)...")
        log_add(log_file = log_file, message = 'Specific error encountered. Retrying with btune = FALSE...')
        # Retry with the tuneStep=c(0.995, 0.9, 0.8, 0.7, 0.6)
        SBayesRC::sbayesrc(
          mafile = paste0(tmp_dir, '/tidy.imp.ma'),
          LDdir = opt$sbayesrc_ldref,
          outPrefix = paste0(tmp_dir, '/sbrc'),
          annot = opt$sbayesrc_annot,
          tuneStep=c(0.995, 0.9, 0.8, 0.7, 0.6),
          log2file = FALSE
        )
      } else {
        # For any other error, rethrow the error and stop the script
      stop("An unexpected error occurred: ", e$message)
    }
  },
  warning = function(w) {
    message("A warning occurred while running SBayesRC::sbayesrc:")
    message(w)
  }
)

###
# Process score file
###

score <- fread(paste0(tmp_dir, '/sbrc.txt'))

# Insert A2 column
ma_imp <- fread(paste0(tmp_dir, '/tidy.imp.ma'))
if(nrow(score) == nrow(ma_imp) & all(score$A1 == ma_imp$A1)){
  score$A2 <- ma_imp$A2
} else {
  log_add(log_file = log_file, message = 'Score file not aligned with sumstats for A2 insertion.')
  system(paste0('cp -r ', tmp_dir, ' ', opt$output_dir,'/temp'))
  stop('Score file not aligned with sumstats for A2 insertion')
}

# Save in plink score format
score <- score[,c('SNP', 'A1', 'A2', 'BETA'), with=F]
names(score) <- c('SNP', 'A1', 'A2', 'SCORE_SBayesRC')

# Flip effects to match reference alleles
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]
score_new <- map_score(ref = ref, score = score)

fwrite(score_new, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

if(file.exists(paste0(opt$output,'.score.gz'))){
  system(paste0('rm ',opt$output,'.score.gz'))
}
system(paste0('gzip ',opt$output,'.score'))

# Record end time of test
if(!is.na(opt$test)){
  test_finish(log_file = log_file, test_start.time = test_start.time)
}

####
# Calculate mean and sd of polygenic scores
####

log_add(log_file = log_file, message = 'Calculating polygenic scores in reference.')

# Calculate scores in the full reference
ref_pgs <- plink_score(pfile = opt$ref_plink_chr, chr = CHROMS, plink2 = opt$plink2, score = paste0(opt$output,'.score.gz'), threads = opt$n_cores)

# Calculate scale within each reference population
pop_data <- read_pop_data(opt$pop_data)

for(pop_i in unique(pop_data$POP)){
  ref_pgs_scale_i <- score_mean_sd(scores = ref_pgs, keep = pop_data[pop_data$POP == pop_i, c('FID','IID'), with=F])
  fwrite(ref_pgs_scale_i, paste0(opt$output, '-', pop_i, '.scale'), row.names = F, quote=F, sep=' ', na='NA')
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()

