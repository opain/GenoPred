#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--config", action="store", default=NULL, type='character',
    help="Pipeline configuration file [required]"),
make_option("--name", action="store", default=NULL, type='character',
    help="Name of target sample [required]"),
make_option("--population", action="store", default=NULL, type='character',
    help="Population in target sample to extract [required]"),
make_option("--plink2", action="store", default='plink2', type='character',
		help="Path PLINK v2 software binary [optional]"),
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
if(is.null(opt$config)){
  stop('--config must be specified.\n')
}
if(is.null(opt$name)){
  stop('--name must be specified.\n')
}
if(is.null(opt$population)){
  stop('--population must be specified.\n')
}

# Read in outdir
outdir <- read_param(config = opt$config, param = 'outdir', return_obj = F)

# Create output directory
opt$output_dir <- paste0(outdir, '/', opt$name, '/pgs/', opt$population)
system(paste0('mkdir -p ',opt$output_dir))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'target_scoring_pipeline.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

# Identify score files to be combined
score_files<-list_score_files(opt$config)

#####
# Combine score files
#####
# They are already sorted, have the same number of variants and have same A1
# Just append without duplicating the SNP, A1, and A2 columns, and insert method and gwas id in column names
scale_files<-list()
for(i in 1:nrow(score_files)){
  if(i == 1){
    system(paste0('zcat ', outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i],".score.gz | sed '1 s/SCORE/",paste0(score_files$method[i],'.',score_files$name[i]),"/g' > ", tmp_dir,'/all_score.txt'))
  } else {
    system(paste0('zcat ', outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i],".score.gz | cut -d' ' --complement -f1-3 - | sed '1 s/SCORE/",paste0(score_files$method[i],'.',score_files$name[i]),"/g' > ", tmp_dir,'/tmp_score.txt'))
    system(paste0("paste -d ' ' ", tmp_dir, '/all_score.txt ', tmp_dir, '/tmp_score.txt > ',tmp_dir, '/temp.txt && mv ', tmp_dir, '/temp.txt ', tmp_dir, '/all_score.txt'))
  }
  # Read in scale file and update Param
  scale_files[[paste0(score_files$method[i],'-',score_files$name[i])]]<-fread(paste0(outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i],'-', opt$population,'.scale'))
  scale_files[[paste0(score_files$method[i],'-',score_files$name[i])]]$Param<-gsub('SCORE', paste0(score_files$method[i],'.',score_files$name[i]), scale_files[[paste0(score_files$method[i],'-',score_files$name[i])]]$Param)
}

# Concatenate scale files
all_scale<-do.call(rbind, scale_files)

#####
# Perform polygenic risk scoring
#####

# Read in target_list
target_list <- read_param(config = opt$config, param = 'target_list', return_obj = T)

# Set params for plink_score
opt$target_plink_chr <- paste0(outdir, '/', opt$name, '/geno/', opt$name, '.ref.chr')
opt$target_keep <- paste0(outdir, '/', opt$name, '/ancestry/keep_files/model_based/', opt$population, '.keep')
refdir <- read_param(config = opt$config, param = 'refdir', return_obj = F)
opt$ref_freq_chr <- paste0(refdir, '/freq_files/', opt$population,'/ref.', opt$population,'.chr')

log_add(log_file = log_file, message = 'Calculating polygenic scores in the target sample.')
scores <-
  plink_score(
    pfile = opt$target_plink_chr,
    chr = CHROMS,
    plink2 = opt$plink2,
    score = paste0(tmp_dir,'/all_score.txt'),
    keep = opt$target_keep,
    frq = opt$ref_freq_chr,
    threads = opt$n_cores
  )

###
# Scale the polygenic scores based on the reference
###

log_add(log_file = log_file, message = 'Scaling target polygenic scores to the reference.')
scores_scaled<-score_scale(score=scores, ref_scale=all_scale)

###
# Write out the target sample scores
###

for(i in 1:nrow(score_files)){
  scores_scaled_i <- scores_scaled[, c('FID','IID', names(scores_scaled)[grepl(paste0(score_files$method[i],'.',score_files$name[i]), names(scores_scaled))]), with=F]
  names(scores_scaled_i) <- gsub(paste0('^', score_files$method[i],'\\.'),'', names(scores_scaled_i))
  dir.create(paste0(outdir, '/', opt$name,'/pgs/', opt$population,'/', score_files$method[i],'/', score_files$name[i]), recursive = T)
  fwrite(scores_scaled_i, paste0(outdir, '/', opt$name,'/pgs/', opt$population,'/', score_files$method[i],'/', score_files$name[i],'/', opt$name,'-', score_files$name[i],'-',opt$population,'.profiles'), sep=' ', na='NA', quote=F)
}

log_add(log_file = log_file, message = paste0('Saved polygenic scores.'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
