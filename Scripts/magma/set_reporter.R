#!/usr/bin/Rscript
# Save start time
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--config", action="store", default=NA, type='character',
              help="config file [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')

# Check required inputs
if(is.null(opt$config)){
  stop('--config must be specified.\n')
}

# Identify outdir from config file
outdir <- read_param(config = opt$config, 'outdir', return_obj = F)

# Read in gwas_list from config file
gwas_list <- read_param(config = opt$config, 'gwas_list', return_obj = T)

set_res<-NULL
for(gwas_i in gwas_list$name){
  if(!file.exists(paste0(outdir,'/reference/gwas_sumstat/',gwas_i,'/magma/sig_indep_sets.txt'))){
    set_res<-rbind(set_res, data.frame(name=gwas_i,
                                       n_sig=0))
  } else {
    
    set_enrich<-read.table(paste0(outdir,'/reference/gwas_sumstat/',gwas_i,'/magma/sig_indep_sets.txt'), header=F)$V1
    
    set_res<-rbind(set_res, data.frame(name=gwas_i,
                                       n_sig=length(set_enrich)))
  }
}

write.table(set_res, paste0(outdir,'/reference/gwas_sumstat/set_reporter.txt'), row.names=F, col.names=T, quote=F)
