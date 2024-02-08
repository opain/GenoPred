#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NULL, type='character',
      help="Path to per chromosome reference PLINK files [required]"),
  make_option("--ldpred2_ref_dir", action="store", default=NULL, type='character',
      help="Path to directory containing LDpred2 reference data [required]"),
  make_option("--ref_keep", action="store", default=NULL, type='character',
      help="Keep file to subset individuals in reference for auto model [optional]"),
  make_option("--pop_data", action="store", default=NULL, type='character',
      help="File containing the population code and location of the keep file [required]"),
  make_option("--plink2", action="store", default='plink2', type='character',
      help="Path PLINK v2 software binary [optional]"),
  make_option("--output", action="store", default=NULL, type='character',
      help="Path for output files [required]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
      help="Number of cores for parallel computing [optional]"),
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
log_header(log_file = log_file, opt = opt, script = 'ldpred2.R', start.time = start.time)

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
sumstats <- read_sumstats(sumstats = opt$sumstats, chr = CHROMS, log_file = log_file, req_cols = c('CHR','SNP','BP','A1','A2','BETA','SE','N','P'))

# Update header for bigsnpr
names(sumstats)<-c('chr','rsid','pos','a1','a0','beta','beta_se','n_eff','p')

# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

# Harmonise with the LDpred2 reference
map<-readRDS(paste0(opt$ldpred2_ref_dir, '/map.rds'))
map<-map[, c('chr', 'pos', 'a0', 'a1', 'af_UKBB', 'ld')]
info_snp <- snp_match(sumstats, map)

#####
# Perform additional suggested QC for LDpred2
#####

# Remove SDss < 0.5 * SDval or SDss > 0.1 + SDval or SDss < 0.1 or SDval < 0.05
sd_val <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))

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

sumstats<-info_snp[!is_bad, ]

log_add(log_file = log_file, message = paste0('Sumstats contains ', nrow(sumstats),' after additional genotype SD check.'))

#######
# Estimate heritability
#######

ldsc <- with(sumstats, snp_ldsc(ld, nrow(map), chi2 = (beta / beta_se)^2, sample_size = n_eff, blocks = NULL))

log_add(log_file = log_file, message = paste0('Estimated SNP-based heritability = ', ldsc[["h2"]]))

if(!is.na(opt$test) & ldsc[["h2"]] < 0.05){
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

beta_inf <- snp_ldpred2_inf(corr, sumstats, ldsc[["h2"]])

log_add(log_file = log_file, message = paste0('Infintesimal model complete at ',as.character(Sys.time())))

#####
# LDpred2-grid
#####

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

####
# LDpred2-auto
####

multi_auto <- snp_ldpred2_auto(corr, sumstats, h2_init = ldsc[["h2"]],
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = opt$n_cores)

# Filter bad chains: `range` should be between 0 and 2
range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))

# Get the final effects using chains that pass this filter
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

log_add(log_file = log_file, message = paste0('Auto model complete at ',as.character(Sys.time())))

####
# Create score file
####

betas <- data.table(SNP=sumstats$rsid, A1=sumstats$a1, A2=sumstats$a0, beta_inf, beta_grid_nosp, beta_grid_sp, beta_auto = beta_auto)

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

fwrite(betas, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

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
ref_pgs <- plink_score(bfile = opt$ref_plink_chr, chr = CHROMS, plink2 = opt$plink2, score = paste0(opt$output,'.score.gz'))

# Calculate scale within each reference population
pop_data <- fread(opt$pop_data)

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
