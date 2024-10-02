#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NULL, type='character',
      help="Path to per chromosome reference PLINK files [required]"),
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
  make_option("--model", action="store", default='grid,auto,inf', type='character',
      help="Specify models to be used by ldpred2 [optional]"),
  make_option("--inference", action="store", default=F, type='character',
      help="Logical indicating whether auto mode should be used for inference [optional]"),
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
log_header(log_file = log_file, opt = opt, script = 'ldpred2.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

# Format the model parameter
opt$model <- unlist(strsplit(opt$model, ','))

# Format the inference parameter
if(!is.logical(opt$inference)){
  opt$inference <- ifelse(opt$inference == 'T', T, F)
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

#######
# Estimate heritability
#######

ldsc <- with(sumstats, snp_ldsc(ld, nrow(map), chi2 = (beta / beta_se)^2, sample_size = n_eff, blocks = NULL))

log_add(log_file = log_file, message = paste0('Estimated SNP-based heritability = ', ldsc[["h2"]]))

if(ldsc[["h2"]] < 0.05){
  ldsc[["h2"]] <- 0.05
}

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
# Run LDpred2
#####

log_add(log_file = log_file, message = 'Running LDpred models.')

# Set seed to ensure reproducibility
set.seed(opt$seed)

#####
# LDpred2-inf
#####

if('inf' %in% opt$model){
  beta_inf <- snp_ldpred2_inf(corr, sumstats, ldsc[["h2"]])

  log_add(log_file = log_file, message = paste0('Infintesimal model complete at ',as.character(Sys.time())))
}

#####
# LDpred2-grid
#####

if('grid' %in% opt$model){

  # Create hyperparameter grid
  h2_seq <- round(ldsc[["h2"]] * c(0.7, 1, 1.4), 4)
  p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
  params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

  beta_grid <- snp_ldpred2_grid(corr, sumstats, params, ncores = opt$n_cores)

  beta_grid_nosp<-data.table(beta_grid[, params$sparse == F])
  names(beta_grid_nosp)<-gsub('-','.',paste0(params$p[params$sparse == F],'_',params$h2[params$sparse == F],'_nosparse'))

  beta_grid_sp<-data.table(beta_grid[,params$sparse == T])
  names(beta_grid_sp)<-gsub('-','.',paste0(params$p[params$sparse == T],'_',params$h2[params$sparse == T],'_sparse'))

  log_add(log_file = log_file, message = paste0('Grid model complete at ',as.character(Sys.time())))

}

####
# LDpred2-auto
####

if('auto' %in% opt$model){
  coef_shrink <- 0.95

  # takes less than 2 min with 4 cores
  multi_auto <- snp_ldpred2_auto(
    corr, sumstats, h2_init = ldsc[["h2"]],
    vec_p_init = seq_log(1e-4, 0.2, length.out = 30), ncores = opt$n_cores,
    allow_jump_sign = FALSE, shrink_corr = coef_shrink)

  # Create convergence plot
  auto <- multi_auto[[1]] # first chain
  png(paste0(opt$output_dir,'/LDpred2_convergence.png'), res=300, unit='px',height=2000, width=2000)
    plot_grid(
      qplot(y = auto$path_p_est) +
        theme_bigstatsr() +
        geom_hline(yintercept = auto$p_est, col = "blue") +
        scale_y_log10() +
        labs(y = "p"),
      qplot(y = auto$path_h2_est) +
        theme_bigstatsr() +
        geom_hline(yintercept = auto$h2_est, col = "blue") +
        labs(y = "h2"),
      ncol = 1, align = "hv"
    )
  dev.off()

  # Filter bad chains: `range` should be between 0 and 2
  range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
  keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))

  # Get the final effects using chains that pass this filter
  beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

  log_add(log_file = log_file, message = paste0('Auto model complete at ',as.character(Sys.time())))
}

####
# Create score file
####

betas <- data.table(SNP=sumstats$rsid, A1=sumstats$a1, A2=sumstats$a0)

if('inf' %in% opt$model){
  betas <- data.table(betas, beta_inf)
}

if('grid' %in% opt$model){
  betas <- data.table(betas, beta_grid_nosp, beta_grid_sp)
}

if('auto' %in% opt$model){
  betas <- data.table(betas, beta_auto = beta_auto)
}

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

####
# Perform inference using LDpred2-auto
####

if(opt$inference){

  log_add(log_file = log_file, message = 'Performing inference...')
  coef_shrink <- 0.95

  multi_auto <- snp_ldpred2_auto(
    corr, sumstats, h2_init = ldsc[["h2"]],
    vec_p_init = seq_log(1e-4, 0.2, length.out = 50), ncores = opt$n_cores,
    burn_in = 500, num_iter = 500, report_step = 20,
    allow_jump_sign = FALSE, shrink_corr = coef_shrink)

  range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
  keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))

  log_add(log_file = log_file, message = '-------------------')
  log_add(log_file = log_file, message = 'Inference results:')

  # Output h2
  all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est, 500))
  auto_h2<-quantile(all_h2, c(0.5, 0.025, 0.975))
  log_add(log_file = log_file, message = paste0('h2 = ', auto_h2[1], '; 95%CI = ', auto_h2[2],' - ', auto_h2[3]))

  # Output p
  all_p <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est, 500))
  auto_p<-quantile(all_p, c(0.5, 0.025, 0.975))
  log_add(log_file = log_file, message = paste0('p = ', auto_p[1], '; 95%CI = ', auto_p[2],' - ', auto_p[3]))

  # Output alpha
  all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500))
  auto_alpha<-quantile(all_alpha, c(0.5, 0.025, 0.975))
  log_add(log_file = log_file, message = paste0('alpha = ', auto_alpha[1], '; 95%CI = ', auto_alpha[2],' - ', auto_alpha[3]))

  # Output PGS-r2
  bsamp <- lapply(multi_auto[keep], function(auto) auto$sample_beta)
  all_r2 <- do.call("cbind", lapply(seq_along(bsamp), function(ic) {
    b1 <- bsamp[[ic]]
    Rb1 <- apply(b1, 2, function(x)
      coef_shrink * bigsparser::sp_prodVec(corr, x) + (1 - coef_shrink) * x)
    b2 <- do.call("cbind", bsamp[-ic])
    b2Rb1 <- as.matrix(Matrix::crossprod(b2, Rb1))
  }))
  auto_r2<-quantile(all_r2, c(0.5, 0.025, 0.975))
  log_add(log_file = log_file, message = paste0('r2 = ', auto_r2[1], '; 95%CI = ', auto_r2[2],' - ', auto_r2[3]))
  log_add(log_file = log_file, message = '-------------------')

  # Output finemapping PIPs
  postp <- rowMeans(sapply(multi_auto[keep], function(auto) auto$postp_est))
  postp <- data.table(SNP=sumstats$rsid, PIP=postp)
  fwrite(postp, paste0(opt$output,'.pip.txt'), col.names=T, sep=' ', quote=F, na='NA')
  log_add(log_file = log_file, message = 'Finemapping results saved.')

}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()
