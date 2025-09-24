#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NA, type='character',
              help="Path to per chromosome reference PLINK files [required]"),
  make_option("--pop_data", action="store", default=NULL, type='character',
              help="File containing the population code and location of the keep file [required]"),
  make_option("--plink2", action="store", default='plink2', type='character',
              help="Path PLINK v2 software binary [optional]"),
  make_option("--output", action="store", default=NULL, type='character',
              help="Path for output files [required]"),
  make_option("--memory", action="store", default=5000, type='numeric',
              help="Memory limit [optional]"),
  make_option("--sumstats", action="store", default=NULL, type='character',
              help="Comma-seperated list of GWAS summary statistics [required]"),
  make_option("--populations", action="store", default=NULL, type='character',
              help="Comma-seperated list of population codes matching GWAS [required]"),
  make_option("--prscsx_path", action="store", default=NULL, type='character',
              help="Path to PRS-CSx executable [required]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
              help="Number of cores for parallel computing [optional]"),
  make_option("--prscsx_ref_path", action="store", default=NULL, type='character',
              help="Comma-seperated list of PRS-CSx reference data [required]"),
  make_option("--test", action="store", default=NA, type='character',
              help="Specify number of SNPs to include [optional]"),
  make_option("--phi_param", action="store", default='auto', type='character',
              help="Path to PRScs reference [optional]"),
  make_option("--seed", action="store", default=1, type='numeric',
              help="Seed number for PRScs [optional]")
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

# Format phi parameters into vector
phi_param<-unlist(strsplit(opt$phi_param,','))
if(any(grepl('auto',phi_param))){
	phi_param<-c(sprintf("%1.00e", as.numeric(phi_param[!grepl('auto',phi_param)])),'auto')
} else {
	phi_param<-sprintf("%1.00e", as.numeric(phi_param))
}

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
if(is.null(opt$prscsx_path)){
  stop('--prscsx_path must be specified.\n')
}
if(is.null(opt$prscsx_ref_path)){
  stop('--prscsx_ref_path must be specified.\n')
}
if(is.null(opt$populations)){
  stop('--populations must be specified.\n')
}

# Create output directory
opt$output_dir <- paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

# Create temp directory
tmp_dir <- tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'prscsx.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

# Split opt$sumstats
sumstats<-unlist(strsplit(opt$sumstats, ','))
log_add(log_file = log_file, message = paste0(length(sumstats), ' sets of GWAS have been provided.'))

gwas_N<-NULL
gwas_CHROMS <- NULL
for(i in 1:length(sumstats)){

  #####
  # Read in sumstats
  #####

  log_add(log_file = log_file, message = 'Reading in GWAS.')

  # Read in, check and format GWAS summary statistics
  gwas <- read_sumstats(sumstats = sumstats[i], chr = CHROMS, log_file = log_file, req_cols = c('CHR','SNP','A1','A2','BETA','SE','N'))
  gwas_CHROMS <- c(gwas_CHROMS, unique(gwas$CHR))
  gwas$CHR<-NULL
  
  # Store average sample size
  gwas_N <- c(gwas_N, round(mean(gwas$N), 0))

  fwrite(gwas, paste0(tmp_dir, '/GWAS_sumstats_',i,'_temp.txt'), sep=' ')

  rm(gwas)
  gc()


}
gwas_CHROMS <- unique(gwas_CHROMS)
print(gwas_CHROMS)
# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

#####
# Process sumstats using PRS-CSx
#####

# Create a temporary reference bim files for PRS-CSx to match to
pvar <- read_pvar(opt$ref_plink_chr, chr = gwas_CHROMS)
pvar$POS<-0
for(i in gwas_CHROMS){
  write.table(pvar[pvar$CHR == i, c('CHR','SNP','POS','BP','A1','A2'), with=F], paste0(tmp_dir,'/ref.chr',i,'.bim'), col.names=F, row.names=F, quote=F)
}

rm(pvar)
gc()

# Make a data.frame listing chromosome and phi combinations
jobs<-NULL
for(i in gwas_CHROMS){
  jobs<-rbind(jobs, data.frame(CHR=i, phi=phi_param))
}

# Run using PRS-CSx auto, and specifying a range of global shrinkage parameters
file.remove(paste0(tmp_dir, '/checker.txt'))
log <- foreach(i = 1:nrow(jobs), .combine = c, .options.multicore = list(preschedule = FALSE)) %dopar% {
  if(!file.exists(paste0(tmp_dir, '/checker.txt'))) {
    # Base command
    command <- paste0(opt$prscsx_path, ' ',
    '--ref_dir=', opt$prscsx_ref_path, '/ ',
    '--bim_prefix=', tmp_dir,'/ref.chr', jobs$CHR[i], ' ',
    '--pop=', opt$populations, ' ',
    '--sst_file=', paste0(paste0(tmp_dir, '/GWAS_sumstats_', 1:length(sumstats),'_temp.txt'), collapse=','),' ',
    '--n_gwas=', paste(gwas_N, collapse=','), ' ',
    '--out_dir=', tmp_dir, '/ ',
    '--out_name=output --chrom=', jobs$CHR[i], ' ',
    '--meta=True --seed=', opt$seed)

    # Add --phi parameter if not 'auto'
    if (jobs$phi[i] != 'auto') {
      command <- paste0(command, ' --phi=', jobs$phi[i])
    }

    # Run command
    log_i <- system(command)

    # Check for an error
    if(log_i != 0){
      write("", paste0(tmp_dir, '/checker.txt'))
    }
  }
}

####
# Combine score files
####

# Read in ref to harmonise score files
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]

score_all<-NULL
for(pop_i in c(unlist(strsplit(opt$populations, ',')), 'META')){
  score_pop<-NULL
  for(phi_i in phi_param){
    score_phi<-NULL
    for(i in gwas_CHROMS){
      score_phi_i<-fread(paste0(tmp_dir,'/output_',pop_i,'_pst_eff_a1_b0.5_phi',phi_i,'_chr',i,'.txt'))
      score_phi<-rbind(score_phi, score_phi_i)
    }
    if(phi_i == phi_param[1]){
      score_phi<-score_phi[,c('V2', 'V4','V5', 'V6'), with=F]
      names(score_phi)<-c('SNP', 'A1', 'A2', paste0('SCORE_',pop_i,'_phi_',phi_i))
    } else {
      score_phi<-score_phi[,'V6', with=F]
      names(score_phi)<-paste0('SCORE_',pop_i,'_phi_',phi_i)
    }
    score_pop<-cbind(score_pop, score_phi)
  }

  # Sort and flip effects to match reference alleles
  score_pop <- map_score(ref = ref, score = score_pop)
  
  if(pop_i == c(unlist(strsplit(opt$populations, ',')), 'META')[1]){
    score_all<-score_pop
  } else {
    score_all<-cbind(score_all, score_pop[, !(names(score_pop) %in% c('SNP','A1','A2')), with=F])
  }
}

fwrite(score_all, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

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
