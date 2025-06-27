#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NULL, type='character',
      help="Path to per chromosome reference PLINK files [required]"),
  make_option("--ref_pcs", action="store", default=NULL, type='character',
      help="Reference PCs for continuous ancestry correction [optional]"),
  make_option("--ldpred2_ref_dir", action="store", default=NULL, type='character',
      help="Path to directory containing LDpred2 reference data [required]"),
  make_option("--pop_data", action="store", default=NULL, type='character',
      help="File containing the population code and location of the keep file [required]"),
  make_option("--plink2", action="store", default='plink2', type='character',
      help="Path PLINK v2 software binary [optional]"),
  make_option("--output", action="store", default=NULL, type='character',
      help="Path for output files [required]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
      help="Number of cores for parallel computing [optional]"),
  make_option("--sample_prev", action="store", default=NULL, type='numeric',
      help="Sampling ratio in GWAS [optional]"),
  make_option("--test", action="store", default=NA, type='character',
      help="Specify number of SNPs to include [optional]"),
  make_option("--binary", action="store", default=F, type='logical',
      help="Specify T if GWAS phenotyp is binary [optional]"),
  make_option("--seed", action="store", default=1, type='numeric',
      help="Set seed to ensure reproducibility  [optional]"),
  make_option("--sumstats", action="store", default=NULL, type='character',
      help="GWAS summary statistics [required]")
)

opt = parse_args(OptionParser(option_list = option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')
library(bigsnpr)
library(ggplot2)

# Check required inputs
if(is.null(opt$ref_plink_chr)){
  stop('--ref_plink_chr must be specified.\n')
}
if(is.null(opt$ldpred2_ref_dir)){
  stop('--ldpred2_ref_dir must be specified.\n')
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
log_header(log_file = log_file, opt = opt, script = 'lassosum2.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

# Format the binary parameter
if(!is.logical(opt$binary)){
  opt$binary <- ifelse(opt$binary == 'T', T, F)
}

if(opt$binary & is.null(opt$sample_prev)){
  stop('--sample_prev must be specified when --binary T.\n')
}

#####
# Read in sumstats
#####

log_add(log_file = log_file, message = 'Reading in GWAS.')

# Read in, check and format GWAS summary statistics
sumstats <- read_sumstats(sumstats = opt$sumstats, chr = CHROMS, log_file = log_file, req_cols = c('CHR','SNP','BP','A1','A2','BETA','SE','N','P'))

# Update header for bigsnpr
names(sumstats)<-c('chr','rsid','pos','a1','a0','beta','beta_se','n_eff','p')

# In binary, update N to be effective N based on opt$sample_prev
if(opt$binary){
  ncas<-sumstats$n_eff*opt$sample_prev
  ncon<-sumstats$n_eff*(1-opt$sample_prev)
  sumstats$n_eff<-4 / (1/ncas + 1/ncon)
  log_add(log_file = log_file, message = paste0('Median effective N = ', median(sumstats$n_eff)))
}

# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

# Harmonise with the LDpred2 reference
map<-readRDS(paste0(opt$ldpred2_ref_dir, '/map.rds'))
names(map)[names(map) == 'af_UKBB']<-'af'
map<-map[, c('chr', 'pos', 'a0', 'a1', 'af', 'ld')]
info_snp <- snp_match(sumstats, map)

#####
# Perform additional suggested QC for LDpred2
#####

# Remove SDss < 0.5 * SDval or SDss > 0.1 + SDval or SDss < 0.1 or SDval < 0.05
sd_val <- with(info_snp, sqrt(2 * af * (1 - af)))

if(opt$binary == F){
  sd_y_est = median(sd_val * info_snp$beta_se * sqrt(info_snp$n_eff))
  sd_ss = with(info_snp, sd_y_est / sqrt(n_eff * beta_se^2))
} else {
  sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))
}

is_bad <-sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05

png(paste0(opt$output_dir,'/LDpred2_sd_qc.png'), res=300, unit='px',height=2000, width=2000)
  plot_obj <- qplot(sd_val, sd_ss, color = is_bad) +
    theme_bigstatsr() +
    coord_equal() +
    scale_color_viridis_d(direction = -1) +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Standard deviations in the validation set",
        y = "Standard deviations derived from the summary statistics",
        color = "Removed?")
  print(plot_obj)
dev.off()

log_add(log_file = log_file, message = paste0('Sumstats contains ', nrow(info_snp[!is_bad, ]),' after additional genotype SD check.'))

sumstats<-info_snp[!is_bad, ]

# If more than half the variants have the wrong SD then the N is probably inaccurate
# Recompute N based on BETA and SE
if(sum(is_bad) > (length(is_bad)*0.5)){
  log_add(log_file = log_file, message = paste0('>50% of variants had a discordant SD.'))
  stop('>50% of variants had a discordant SD. Check the sample size in the sumstats.')
}

#####
# Prepare LD reference data
#####

log_add(log_file = log_file, message = 'Creating genome-wide sparse matrix.')

# Create genome-wide sparse LD matrix
for (chr in CHROMS) {
  ## indices in 'sumstats'
  ind.chr <- which(sumstats$chr == chr)
  ## indices in 'map'
  ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map$chr == chr))
  
  corr0 <- readRDS(paste0(opt$ldpred2_ref_dir, '/LD_with_blocks_chr', chr, '.rds'))[ind.chr3, ind.chr3]
  
  if (chr == CHROMS[1]) {
    corr <- as_SFBM(corr0, paste0(tmp_dir, '/LD_GW_sparse'), compact = TRUE)
  } else {
    corr$add_columns(corr0, nrow(corr))
  }
}

#####
# Run lassosum2
#####

log_add(log_file = log_file, message = 'Running lassosum2.')

# Set seed to ensure reproducibility
set.seed(opt$seed)

# lassosum2
beta_df <- snp_lassosum2(
  corr,
  sumstats,
  ncores = opt$n_cores
)

####
# Create score file
####

# Convert matrix to data.table with column names
beta_dt <- as.data.table(beta_df)
grid <- attr(beta_df, "grid_param")
new_names <- paste0("s", grid$delta, "_lambda", grid$lambda)
setnames(beta_dt, new_names)

betas <- data.table(SNP=sumstats$rsid, A1=sumstats$a1, A2=sumstats$a0, beta_dt)

rem<-NULL
for(i in 4:length(names(betas))){
  if(is.infinite(sum(betas[[names(betas)[i]]])) | is.na(sum(betas[[names(betas)[i]]]))){
    log_add(log_file = log_file, message = paste0('Skipping ',names(betas)[i],' due to presence of non-finite values.'))
    rem<-c(rem,i)
  }
}

if(is.null(rem) == F){
  betas<-betas[, -rem, with=F]
}

names(betas)[-1:-3] <- paste0('SCORE_', names(betas)[-1:-3])

# Flip effects to match reference alleles
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]
score_new <- map_score(ref = ref, score = betas)

# Reduce number of significant figures to save space
score_new[, (4:ncol(score_new)) := lapply(.SD, signif, digits = 7), .SDcols = 4:ncol(score_new)]

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
