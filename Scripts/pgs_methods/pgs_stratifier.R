#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--config", action="store", default=NULL, type='character',
    help="Pipeline configuration file [required]"),
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
library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

# Check required inputs
if(is.null(opt$config)){
  stop('--config must be specified.\n')
}

# Read in outdir
outdir <- read_param(config = opt$config, param = 'outdir', return_obj = F)

# Read in refdir
refdir <- read_param(config = opt$config, param = 'refdir', return_obj = F)

# Read in resdir
resdir <- read_param(config = opt$config, param = 'resdir', return_obj = F)

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(outdir, '/reference/pgs_score_files/stratifier_', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), '.log')
log_header(log_file = log_file, opt = opt, script = 'target_scoring_partitioned_pipeline.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

# Identify score files to be combined
score_files<-list_score_files(opt$config)

# Restrict to single source PGS
score_files <- score_files[!(score_files$method %in% pgs_group_methods) & !grepl('tlprs|leopard', score_files$method),]

# Check which score files need to be partitioned score files or target genetic data are newer than target pgs
if(!is.null(score_files)){
  set_reporter_file <- paste0(outdir, '/reference/gwas_sumstat/set_reporter.txt')
  set_reporter<-fread(set_reporter_file)
  set_reporter_file_time <- file.info(set_reporter_file)$mtime

  # Remove score files for gwas that have no significant sets
  score_files<-score_files[score_files$name %in% set_reporter$name[set_reporter$n_sig > 0],]
  
  score_files_to_do <- data.table()
  for(i in 1:nrow(score_files)){
    score_i <- paste0(outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i], '.score.gz')
    score_partitioned_i <- paste0(outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i], '.partitioned.score.gz')
    if(!file.exists(score_partitioned_i)){
      score_files_to_do <- rbind(score_files_to_do, score_files[i,])
    } else {
      score_i_time <- file.info(score_i)$mtime
      score_partitioned_i_time <- file.info(score_partitioned_i)$mtime
      if (score_i_time > pgs_i_time | set_reporter_file_time > score_partitioned_i_time) {
        score_files_to_do <- rbind(score_files_to_do, score_files[i,])
        system(paste0('rm ', pgs_i))
      }
    }
  }
  log_add(log_file = log_file, message = paste0('After checking timestamps, ', nrow(score_files_to_do), '/', nrow(score_files), ' score files will be partitioned.'))
  score_files <- score_files_to_do
}

if(is.null(score_files) || nrow(score_files) == 0){
  log_add(log_file = log_file, message = paste0('No score files to be partitioned.'))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Analysis finished at',as.character(end.time),'\n')
  cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  quit(save = "no", status = 0)
}

# Stratify score files (pseudo only)
log <- foreach(i = 1:nrow(score_files), .combine = c, .options.multicore = list(preschedule = FALSE)) %dopar% {
  # Read in snplists
  set_enrich<-read.table(paste0(outdir,'/reference/gwas_sumstat/',score_files$name[i],'/magma/sig_indep_sets.txt'), header=F)$V1

  param <- find_pseudo(
    config = opt$config,
    gwas = score_files$name[i],
    pgs_method = score_files$method[i]
  )
    
  score_header <- fread(
    paste0(outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-', score_files$name[i],".score.gz"), 
    nrows = 1)
  score_cols <- which(names(score_header) %in% c('SNP', 'A1', 'A2', paste0('SCORE_', param)))
  
  # Create stratified score files
  score_i <- fread(cmd = paste0(
    'zcat ', outdir, '/reference/pgs_score_files/', score_files$method[i], '/', score_files$name[i], '/ref-', score_files$name[i], ".score.gz | ",
    "cut -d' ' -f ", paste(score_cols, collapse =','), " - ")) # Keep pseudo score
  
  for(k in 1:length(set_enrich)){
    snplist_k <- fread(paste0(outdir,'/reference/gwas_sumstat/',score_files$name[i],'/magma/snplists/',set_enrich[k],'.snplist'), header=F)$V1
    score_i[[paste0(names(score_i)[4], '.set_', k)]] <- score_i[[4]]
    score_i[[paste0(names(score_i)[4], '.set_', k)]][!(score_i$SNP %in% snplist_k)] <- 0
  }
  
  # Remove unstratified PGS
  score_i<-score_i[,-4]
  
  # Save score file
  file_name <- paste0(outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-', score_files$name[i],".stratified")
  fwrite(score_i, 
         paste0(file_name, '.score'), 
         col.names=T, sep=' ', quote=F)
  
  if(file.exists(paste0(file_name, ".score.gz"))){
    system(paste0('rm ', file_name, ".score.gz"))
  }
  system(paste0('gzip ', file_name, '.score'))
  
  # Calculate scores in the full reference
  # Subset to variants with non-zero effect
  extract_snplist <- score_i$SNP[rowSums(abs(score_i[,-1:-3])) != 0]
  ref_pgs <-
    plink_score(
      pfile = paste0(refdir, '/ref.chr'),
      chr = CHROMS,
      plink2 = opt$plink2,
      extract = extract_snplist,
      score = paste0(file_name, '.score.gz'),
      threads = opt$n_cores
    )
  
  # Derive trans-ancestry PGS models and estimate PGS residual scale
  model_trans_pgs(
    scores = ref_pgs,
    pcs = paste0(resdir, '/data/ref/pc_score_files/TRANS/ref-TRANS-pcs.profiles'),
    output = file_name
  )
  
  # Calculate scale within each reference population
  pop_data <- read_pop_data(paste0(refdir, '/ref.pop.txt'))
  
  for(pop_i in unique(pop_data$POP)){
    ref_pgs_scale_i <- score_mean_sd(scores = ref_pgs, keep = pop_data[pop_data$POP == pop_i, c('FID','IID'), with=F])
    fwrite(ref_pgs_scale_i, paste0(file_name, '-', pop_i, '.scale'), row.names = F, quote=F, sep=' ', na='NA')
  }
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
