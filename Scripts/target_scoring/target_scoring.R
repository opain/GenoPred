#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--target_plink_chr", action="store", default=NULL, type='character',
		help="Path to per chromosome target PLINK2 files [required]"),
make_option("--target_keep", action="store", default=NULL, type='character',
		help="Path to keep file for target [optional]"),
make_option("--ref_score", action="store", default=NULL, type='character',
		help="Path to reference scoring files [required]"),
make_option("--ref_freq_chr", action="store", default=NULL, type='character',
		help="Path to per chromosome reference PLINK2 .afreq files [required]"),
make_option("--plink2", action="store", default='plink2', type='character',
		help="Path PLINK v2 software binary [optional]"),
make_option("--output", action="store", default=NULL, type='character',
		help="Path for output files [required]"),
make_option("--ref_scale", action="store", default=NULL, type='character',
		help="Path reference scale file [required]"),
make_option("--pheno_name", action="store", default=NULL, type='character',
    help="Name of phenotype to be added to column names. Default is SCORE. [optional]"),
make_option("--n_cores", action="store", default=1, type='numeric',
    help="Number of cores to use [optional]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')

# Check required inputs
if(is.null(opt$target_plink_chr)){
  stop('--target_plink_chr must be specified.\n')
}
if(is.null(opt$ref_score)){
  stop('--ref_score must be specified.\n')
}
if(is.null(opt$ref_freq_chr)){
  stop('--ref_freq_chr must be specified.\n')
}
if(is.null(opt$ref_score)){
  stop('--ref_score must be specified.\n')
}
if(is.null(opt$output)){
  stop('--output must be specified.\n')
}
if(is.null(opt$ref_scale)){
  stop('--ref_scale must be specified.\n')
}

# Create output directory
opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'target_scoring.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

#####
# Perform polygenic risk scoring
#####

log_add(log_file = log_file, message = 'Calculating polygenic scores in the target sample.')
scores<-plink_score(bfile = opt$target_plink_chr, chr = CHROMS, plink2 = opt$plink2, score = opt$ref_score, keep = opt$target_keep, frq = opt$ref_freq_chr, threads = opt$n_cores)

###
# Scale the polygenic scores based on the reference
###

log_add(log_file = log_file, message = 'Scaling target polygenic scores to the reference.')
ref_scale<-fread(opt$ref_scale)
scores_scaled<-score_scale(score=scores, ref_scale=ref_scale)

###
# Write out the target sample scores
###

if(!is.null(opt$pheno_name)){
	names(scores_scaled)<-gsub('SCORE', opt$pheno_name, names(scores_scaled))
}

fwrite(scores_scaled, paste0(opt$output,'.profiles'), sep=' ', na='NA', quote=F)
log_add(log_file = log_file, message = paste0('Saved polygenic scores to: ',opt$output,'.profiles.'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
