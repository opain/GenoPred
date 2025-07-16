#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NULL, type='character',
      help="Path to per chromosome reference PLINK files [required]"),
  make_option("--ref_pcs", action="store", default=NULL, type='character',
      help="Reference PCs for continuous ancestry correction [optional]"),
  make_option("--ref_freq_chr", action="store", default=NULL, type='character',
      help="Path to per chromosome reference PLINK .frq files [required]"),
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
  make_option("--impute_N", action="store", default=T, type='logical',
      help="Logical indicating whether per variant N should imputed based on SE. [optional]"),
  make_option("--P_max", action="store", default=NULL, type='numeric',
      help="P-value threshold for filter variants [optional]"),
  make_option("--robust", action="store", default=F, type='logical',
      help="Force robust GCTB parameterisation [optional]"),
  make_option("--test", action="store", default=NA, type='character',
      help="Specify number of SNPs to include [optional]"),
  make_option("--pop_prev", action="store", default=NULL, type='numeric',
      help="Population prevelance (if binary) [optional]"),
  make_option("--sample_prev", action="store", default=NULL, type='numeric',
      help="Sampling ratio in GWAS [optional]"),
  make_option("--ld_matrix_chr", action="store", default=NULL, type='character',
      help="Path to per chromosome shrunk sparse LD matrix from GCTB [required]")
)

opt = parse_args(OptionParser(option_list = option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')
library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

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
if(is.null(opt$ld_matrix_chr)){
  stop('--ld_matrix_chr must be specified.\n')
}
if(is.na(as.numeric(opt$sample_prev))){
  opt$sample_prev<-NULL
}
if(is.na(as.numeric(opt$pop_prev))){
  opt$pop_prev<-NULL
}
if(any(!is.null(c(opt$sample_prev, opt$pop_prev))) & any(is.null(c(opt$sample_prev, opt$pop_prev)))){
  stop('If either sample_prev or pop_prev are specified, both must be specified.')
}

# Create output directory
opt$output_dir <- paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

# Create temp directory
tmp_dir <- tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'sbayesr.R', start.time = start.time)

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
gwas <- read_sumstats(sumstats = opt$sumstats, chr = CHROMS, log_file = log_file, req_cols = c('CHR','SNP','A1','A2','FREQ','BETA','SE','P','N'))
GWAS_CHROMS<-unique(gwas$CHR)
gwas$CHR<-NULL
###
# Change to COJO format
###

gwas <- gwas[, c('SNP','A1','A2','FREQ','BETA','SE','P','N'), with=F]
names(gwas) <- c('SNP','A1','A2','freq','b','se','p','N')

# Check whether per variant sample size is available
if(length(unique(gwas$N)) == 1){
  per_var_N <- F
  log_add(log_file = log_file, message = 'Per variant N is not present.')

  if(opt$impute_N == T){
    log_add(log_file = log_file, message = 'Per variant N will be imputed.')
  }
} else {
  per_var_N <- T
  log_add(log_file = log_file, message = 'Per variant N is present.')
}

# Set maximum p-value threshold
if(!is.null(opt$P_max)){
  gwas <- gwas[gwas$p <= opt$P_max,]
  log_add(log_file = log_file, message = paste0('After p-value threshold of <= ',opt$P_max,', ', nrow(gwas), ' variants remain.'))
}

# Write out cojo format sumstats
fwrite(gwas, paste0(tmp_dir,'/GWAS_sumstats_COJO.txt'), sep=' ', na = "NA", quote=F)

rm(gwas)
gc()

# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

#####
# Run GCTB SBayesR
#####

log_add(log_file = log_file, message = 'Running SBayesR analysis.')

sbayesr_opt <- NULL
if(opt$robust){
  sbayesr_opt <- paste0(sbayesr_opt, '--robust ')
}
if(per_var_N == F & opt$impute_N == T){
  sbayesr_opt <- paste0(sbayesr_opt, '--impute-n ')
}

error<-foreach(i = GWAS_CHROMS, .combine = rbind, .options.multicore = list(preschedule = FALSE)) %dopar% {
  log <- system(paste0(opt$gctb, ' --sbayes R --ldm ', opt$ld_matrix_chr, i, '.ldm.sparse --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1 --gwas-summary ', tmp_dir, '/GWAS_sumstats_COJO.txt --chain-length 10000 ', sbayesr_opt, '--exclude-mhc --burn-in 2000 --out-freq 1000 --out ', tmp_dir, '/GWAS_sumstats_SBayesR.chr', i),  intern = T)

  # Check whether the analysis converged
  if(any(grepl("Analysis finished", log))){
    if(any(grepl("MCMC cycles completed", log))){
	    data.frame(chr=i, Log='Analysis converged')
    } else {
      print(log)
      data.frame(chr=i, Log='Analysis did not converge')
    }
  } else {
    print(log)
	  data.frame(chr=i, Log='Error')
  }
}

# Report an error if SBayesR didn't converge for all chromosomes
if(sum(grepl('Error', error$Log) == T) > 1){
  log_add(log_file = log_file, message = paste0('An error occurred for ', sum(grepl('Error', error$Log) == T), ' chromosomes. Retry requesting more memory or run interactively to debug.'))
  sink(file = log_file, append = T)
    print(error)
  sink()
  q()
  n
}

# Combine per chromosome snpRes files
snpRes<-NULL
for(i in GWAS_CHROMS){
  snpRes <- rbind(snpRes, fread(paste0(tmp_dir, '/GWAS_sumstats_SBayesR.chr', i, '.snpRes')))
}

# Save in plink score format
snpRes <- snpRes[,c('Name', 'A1', 'A2', 'A1Effect'), with=F]
names(snpRes) <- c('SNP', 'A1', 'A2', 'SCORE_SBayesR')

# Flip effects to match reference alleles
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]
score_new <- map_score(ref = ref, score = snpRes)

fwrite(score_new, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

if(file.exists(paste0(opt$output,'.score.gz'))){
  system(paste0('rm ',opt$output,'.score.gz'))
}
system(paste0('gzip ',opt$output,'.score'))

# Record end time of test
if(!is.na(opt$test)){
  test_finish(log_file = log_file, test_start.time = test_start.time)
}

# Combine per chromosome parRes files
parRes_mcmc <- list()
for(i in GWAS_CHROMS){
	parRes_mcmc[[i]] <- fread(paste0(tmp_dir, '/GWAS_sumstats_SBayesR.chr', i, '.mcmcsamples.Par'))
}

parRes <- NULL
for(par in names(parRes_mcmc[[i]])){
	parRes_mcmc_par <- NULL
	for(i in GWAS_CHROMS){
		parRes_mcmc_par <- cbind(parRes_mcmc_par, parRes_mcmc[[i]][[par]])
	}

	parRes_mcmc_par_sum <- rowSums(parRes_mcmc_par)

	parRes_par <- data.frame( Par = par,
                            Mean = mean(parRes_mcmc_par_sum),
                            SD = sd(parRes_mcmc_par_sum))

	parRes <- rbind(parRes, parRes_par)
}

write.table(parRes, paste0(opt$output_dir,'/GWAS_sumstats_SBayesR.GW.parRes'), col.names=T, row.names=F, quote=F)
log_add(log_file = log_file, message = paste0('SNP-heritability estimate is ',parRes[parRes$Par == 'hsq', names(parRes) == 'Mean']," (SD=",parRes[parRes$Par == 'hsq', names(parRes) == 'SD'],")."))

write.table(parRes[parRes$Par == 'hsq', names(parRes) == 'Mean'], paste0(opt$output,'.hsq_obs'), col.names=F, row.names = F, quote = F)

# Convert the SNP-heritability to the liability scale
if(!is.null(opt$pop_prev) && !is.null(opt$sample_prev)) {
  h2_liab <- h2l_R2(k = opt$pop_prev, r2 = parRes[parRes$Par == 'hsq', names(parRes) == 'Mean'], p = opt$sample_prev)
  h2_liab_se <- h2l_R2(k = opt$pop_prev, r2 = parRes[parRes$Par == 'hsq', names(parRes) == 'SD'], p = opt$sample_prev)
  log_add(log_file = log_file, message = paste0('SNP-heritability estimate on the liability scale = ', round(h2_liab, 4), " (", round(h2_liab_se, 4), ")."))
  
  write.table(h2_liab, paste0(opt$output,'.hsq_liab'), col.names=F, row.names = F, quote = F)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()

