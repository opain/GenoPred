#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--config", action="store", default=NULL, type='character',
              help="Pipeline configuration file [required]"),
  make_option("--test", action="store", default=NA, type='character',
              help="Specify number of SNPs to include [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')

# Read in outdir
outdir <- read_param(config = opt$config, param = 'outdir', return_obj = F)

# Create output directory
opt$output <- paste0(outdir, '/reference/ref')
system(paste0('mkdir -p ',opt$output))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output,'/finding_intersecting_variants.log')
log_header(log_file = log_file, opt = opt, script = 'find_intersecting_variants.R', start.time = start.time)

# Check required inputs
if(is.null(opt$config)){
  stop('--config must be specified.\n')
}

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

########
# Read in pvar files for all target datasets
########

target_list <- read_param(config = opt$config, param = 'target_list', return_obj = T)

log_add(log_file = log_file, message = paste0('Reading pvar files from ', nrow(target_list),' target datasets.'))

pvar_list<-list()
for(target in target_list$name){
  pvar_list[[target]]<-read_pvar(paste0(outdir, '/', target,'/geno/',target,'.ref.chr'), chr = CHROMS)$SNP
}
intersecting_variants <- Reduce(intersect, pvar_list)

log_add(log_file = log_file, message = paste0(length(intersecting_variants),' variants in common across target datasets.'))

#######
# Subset the reference data to variants present in all target datasets
#######

log_add(log_file = log_file, message = paste0('Creating subset of reference data restricted to variants in target datasets...'))

refdir <- read_param(config = opt$config, param = 'refdir', return_obj = F)

plink_subset(
  plink2 = 'plink2',
  pfile = paste0(refdir, '/ref.chr'),
  chr = CHROMS,
  extract = intersecting_variants,
  out = paste0(opt$output, '/ref.chr')
)

for(chr in CHROMS){
  rds<-readRDS(paste0(refdir, '/ref.chr', chr,'.rds'))
  rds<-rds[rds$SNP %in% intersecting_variants,]
  saveRDS(rds, paste0(opt$output, '/ref.chr', chr, '.rds'))
}

log_add(log_file = log_file, message = paste0('Reference data subset saved to ', opt$output, '/ref.chr'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()
