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
library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

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
opt$output <- paste0(outdir, '/', opt$name, '/pgs/', opt$population)
system(paste0('mkdir -p ',opt$output))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output, '_partitioned_', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), '.log')
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

# Check whether score files or target genetic data are newer than target pgs
if(!is.null(score_files)){
  ancestry_reporter_file<-paste0(outdir, '/reference/target_checks/', opt$name, '/ancestry_reporter.done')
  ancestry_reporter_file_time <- file.info(ancestry_reporter_file)$mtime
  
  set_reporter_file <- paste0(outdir, '/reference/gwas_sumstat/set_reporter.txt')
  set_reporter<-fread(set_reporter_file)

  # Remove score files for gwas that have no significant sets
  score_files<-score_files[score_files$name %in% set_reporter$name[set_reporter$n_sig > 0],]
  
  score_files_to_do <- data.table()
  for(i in 1:nrow(score_files)){
    pgs_i <- paste0(outdir, '/', opt$name,'/pgs/', opt$population,'/', score_files$method[i],'/', score_files$name[i],'/', opt$name,'-', score_files$name[i],'-',opt$population,'.partitioned.profiles')
    score_i <- paste0(outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i], '.stratified.score.gz')
    if(!file.exists(pgs_i)){
      score_files_to_do <- rbind(score_files_to_do, score_files[i,])
    } else {
      score_i_time <- file.info(score_i)$mtime
      pgs_i_time <- file.info(pgs_i)$mtime
      if (score_i_time > pgs_i_time | ancestry_reporter_file_time > pgs_i_time) {
        score_files_to_do <- rbind(score_files_to_do, score_files[i,])
        system(paste0('rm ', pgs_i))
      }
    }
  }
  log_add(log_file = log_file, message = paste0('After checking timestamps, ', nrow(score_files_to_do), '/', nrow(score_files), ' score files will be used for target scoring.'))
  score_files <- score_files_to_do
}

if(is.null(score_files) || nrow(score_files) == 0){
  log_add(log_file = log_file, message = paste0('No score files to be used for target scoring.'))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  sink(file = log_file, append = T)
  cat('Analysis finished at',as.character(end.time),'\n')
  cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  quit(save = "no", status = 0)
}

# Set params for plink_score
opt$target_plink_chr <- paste0(outdir, '/', opt$name, '/geno/', opt$name, '.ref.chr')
if(opt$population == 'TRANS'){
  opt$target_keep<-NULL
} else {
  opt$target_keep <- paste0(outdir, '/', opt$name, '/ancestry/keep_files/model_based/', opt$population, '.keep')
}
refdir <- read_param(config = opt$config, param = 'refdir', return_obj = F)
opt$ref_freq_chr <- paste0(refdir, '/freq_files/', opt$population,'/ref.', opt$population,'.chr')

# Read in reference SNP data
ref <- read_pvar(paste0(refdir, '/ref.chr'), chr = CHROMS)[, c('CHR','SNP','A1','A2'), with=F]

# Identify SNPs within sets
set_snps<-NULL
for(i in unique(score_files$name)){
  set_enrich<-read.table(paste0(outdir,'/reference/gwas_sumstat/',i,'/magma/sig_indep_sets.txt'), header=F)$V1
  for(j in set_enrich){
    set_snps <- c(set_snps, fread(paste0(outdir,'/reference/gwas_sumstat/',i,'/magma/snplists/',j,'.snplist'), header=F)$V1)
  }
}
set_snps <- unique(set_snps)
ref$extract <- ifelse(ref$SNP %in% set_snps, T, F)

# We will process score files and perform target scoring for one chromosome for efficiency
for(chr_i in CHROMS){
  log_add(log_file = log_file, message = '########################')
  log_add(log_file = log_file, message = paste0('Processing chromosome ', chr_i,':'))

  #####
  # Combine score files
  #####
  # Only retain pseudo score
  
  # Create row number index to subset score files by chromosome
  row_index <- format(which(ref$CHR == chr_i & ref$extract == T) + 1, scientific = FALSE)
  write.table(row_index, paste0(tmp_dir,'/row_index.txt'), row.names=F, quote=F, col.names = F)
  ref_subset <- ref[as.numeric(row_index),]
  
  # Create file containing SNP, A1, and A2 information for each chromosome
  fwrite(ref_subset[, c('SNP','A1','A2'), with=F], paste0(tmp_dir,'/map.txt'), row.names=F, quote=F, sep=' ')

  # Extract process score files for each name (gwas/score) in parallel
  foreach(i = 1:nrow(score_files), .combine = c, .options.multicore = list(preschedule = FALSE)) %dopar% {
    system(paste0(
      'zcat ', outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i],".stratified.score.gz | ",
      "cut -d' ' --complement -f1-3 | ",  # Keep relevant columns, remove first 3
      'awk \'NR==FNR {rows[$1]; next} FNR==1 || FNR in rows\' ', paste0(tmp_dir,'/row_index.txt'), ' - | ',  # Corrected to retain the header and process indexed rows
      "sed '1 s/SCORE_/", paste0('score_file_', i,'.'), "/g' > ",  # Replace SCORE in the header
      tmp_dir, '/tmp_score.', paste0(score_files$method[i], '.', score_files$name[i]), '.txt'
    ))
  }

  # Paste files together in batches
  # Set number of batches according to the number of score files to combine
  num_batches <- max(c(1, min(c(opt$n_cores, floor(nrow(score_files) / 2)))))
  tmp_score_files <- paste0(tmp_dir,'/tmp_score.',score_files$method,'.',score_files$name,'.txt')
  set.seed(1)
  batches <- split(sample(tmp_score_files), rep(1:num_batches, length.out = length(tmp_score_files)))
  log_add(log_file = log_file, message = paste0('Aggregating score files in ', num_batches,' batches.'))
  foreach(i = 1:length(batches), .combine = c, .options.multicore = list(preschedule = FALSE)) %dopar% {
    system(paste0("paste -d ' ' ", paste(batches[[i]], collapse = " "),' > ',tmp_dir,'/tmp_batch_',i))
    system(paste0('rm ', paste(batches[[i]], collapse = " ")))
  }

  # Paste batches together
  log_add(log_file = log_file, message = paste0('Aggregating batched score files.'))
  tmp_batch_files <- paste0(tmp_dir,'/tmp_batch_',1:length(batches))
  system(paste0("paste -d ' ' ", tmp_dir,'/map.txt ', paste(tmp_batch_files, collapse = " "), ' > ', tmp_dir, '/all_score.txt'))
  system(paste0('rm ', paste(tmp_batch_files, collapse = " ")))

  # Perform polygenic risk scoring
  scores_i <-
    plink_score(
      pfile = opt$target_plink_chr,
      chr = chr_i,
      plink2 = opt$plink2,
      score = paste0(tmp_dir,'/all_score.txt'),
      keep = opt$target_keep,
      frq = opt$ref_freq_chr,
      threads = opt$n_cores
    )

  # Sum scores across chromosomes
  if(chr_i == CHROMS[1]){
    scores_ids <- scores_i[, 1:2, with = F]
    current_scores <- as.matrix(scores_i[, -1:-2, with = FALSE])
    scores <- current_scores
  } else {
    current_scores <- as.matrix(scores_i[, -1:-2, with = FALSE])
    scores <- scores + current_scores
  }

  system(paste0('rm ', tmp_dir, '/all_score.txt'))
  system(paste0('rm ', tmp_dir, '/row_index.txt'))
  system(paste0('rm ', tmp_dir, '/map.txt'))
}

# Combine score with IDs
scores<-data.table(scores_ids,
                   scores)

###
# Scale the polygenic scores based on the reference
###

if(opt$population == 'TRANS'){
  log_add(log_file = log_file, message = paste0('Reading in ancestry adjustment models.'))

  models<-list()
  for(i in 1:nrow(score_files)){
    models[[paste0('score_file_', i)]]<-readRDS(paste0(outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i],'.stratified-TRANS.model.rds'))
    names(models[[paste0('score_file_', i)]])<-gsub('SCORE_', paste0('score_file_', i, '.'), names(models[[paste0('score_file_', i)]]))
  }
  
  models <- do.call(c, unname(models))
  
  # Read in target projected PCs
  target_pcs<-fread(paste0(outdir,'/',opt$name,'/pcs/projected/TRANS/',opt$name,'-TRANS.profiles'))
  log_add(log_file = log_file, message = paste0('Reading in target reference-projected PCs.'))
  
  # Adjust scores
  log_add(log_file = log_file, message = 'Adjusting target PGS for ancestry.')
  scores <- score_adjust(score = scores, pcs = target_pcs, ref_model = models)
} else {
  # Read in scale file and update Param
  log_add(log_file = log_file, message = paste0('Reading in scale files.'))
  scale_files<-list()
  for(i in 1:nrow(score_files)){
    scale_files[[paste0('score_file_', i)]]<-fread(paste0(outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i],'.stratified-', opt$population,'.scale'))
    scale_files[[paste0('score_file_', i)]]$Param<-gsub('SCORE_', paste0('score_file_', i, '.'), scale_files[[paste0('score_file_', i)]]$Param)
  }
  
  # Concatenate scale files
  all_scale<-do.call(rbind, scale_files)
  
  # Scale scores
  log_add(log_file = log_file, message = 'Scaling target polygenic scores to the reference.')
  scores<-score_scale(score=scores, ref_scale=all_scale)
}

###
# Write out the target sample scores
###

for(i in 1:nrow(score_files)){
  scores_i <- scores[, c('FID','IID', names(scores)[grepl(paste0('^score_file_', i, '\\.'), names(scores))]), with=F]
  names(scores_i) <- gsub(paste0('^score_file_', i, '\\.'), paste0(score_files$name[i], '_'), names(scores_i))
  dir.create(paste0(outdir, '/', opt$name,'/pgs/', opt$population,'/', score_files$method[i],'/', score_files$name[i]), recursive = T)
  fwrite(scores_i, paste0(outdir, '/', opt$name,'/pgs/', opt$population,'/', score_files$method[i],'/', score_files$name[i],'/', opt$name,'-', score_files$name[i],'-',opt$population,'.partitioned.profiles'), sep=' ', na='NA', quote=F)
}

log_add(log_file = log_file, message = paste0('Saved polygenic scores.'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
