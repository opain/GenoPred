#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--config", action="store", default=NULL, type='character',
    help="Pipeline configuration file [required]"),
make_option("--continuous", action="store", default=T, type='logical',
    help="Logical indicating whether or not continuous correction for ancestry is required [optional]"),
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
library(glm2)
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

# Create output directory
opt$output <- paste0(outdir, '/reference/pgs_score_files/ref_scoring')

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output, '_', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), '.log')
log_header(log_file = log_file, opt = opt, script = 'ref_scoring_pipeline.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

# Identify score files to be combined
score_files<-list_score_files(opt$config)
refdir <- read_param(config = opt$config, param = 'refdir', return_obj = F)
ref <- read_pvar(paste0(refdir, '/ref.chr'), chr = CHROMS)[, c('CHR','SNP','A1','A2'), with=F]
reference_fingerprint <- reference_variant_fingerprint(ref)
if(!is.null(score_files)) for(i in seq_len(nrow(score_files))) {
  validate_score_reference(paste0(outdir, '/reference/pgs_score_files/', score_files$method[i], '/', score_files$name[i], '/ref-', score_files$name[i], '.score.gz'), ref, reference_fingerprint)
}

# Check whether score files or target genetic data are newer than target pgs
if(!is.null(score_files)){
  ref_pcs_files <- character()
  if(opt$continuous){
    ref_pcs_files <- paste0(outdir, '/reference/pc_score_files/TRANS/ref-TRANS-pcs', c('.EUR.scale', '.profiles'))
  }
  
  score_files_to_do <- data.table()
  for(i in 1:nrow(score_files)){
    pgs_i <- paste0(outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i], '-EUR.profiles')
    score_i <- paste0(outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i], '.score.gz')
    if(reference_score_needs_update(score_i, sub('-EUR.profiles$', '', pgs_i), opt$continuous, ref_pcs_files)){
      score_files_to_do <- rbind(score_files_to_do, score_files[i,])
    }
  }
  log_add(log_file = log_file, message = paste0('After checking timestamps, ', nrow(score_files_to_do), '/', nrow(score_files), ' score files will be used for reference scoring.'))
  score_files <- score_files_to_do
}

if(is.null(score_files) || nrow(score_files) == 0){
  log_add(log_file = log_file, message = paste0('No score files to be used for reference scoring.'))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  sink(file = log_file, append = T)
  cat('Analysis finished at',as.character(end.time),'\n')
  cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  quit(save = "no", status = 0)
}

# Set params for plink_score
refdir <- read_param(config = opt$config, param = 'refdir', return_obj = F)
opt$ref_plink_chr <- paste0(refdir, '/ref.chr')

# Read in reference SNP data
# Reference SNP data was loaded above to validate even already-scored artifacts.

ref_map <- paste0(tmp_dir, '/score_reference.tsv')
fwrite(ref[, c('CHR','SNP','A1','A2'), with=FALSE], ref_map, sep='\t', col.names=FALSE, quote=FALSE)
for(i in seq_len(nrow(score_files))) {
  score_path <- paste0(outdir, '/reference/pgs_score_files/', score_files$method[i], '/', score_files$name[i], '/ref-', score_files$name[i], '.score.gz')
  prepare_score_shards(score_path, ref_map, paste0(tmp_dir, '/score_shard.', i), i)
}

# We will process score files and perform target scoring for one chromosome for efficiency
for(chr_i in CHROMS){
  log_add(log_file = log_file, message = '########################')
  log_add(log_file = log_file, message = paste0('Processing chromosome ', chr_i,':'))

  #####
  # Combine score files
  #####
  # Create file containing SNP, A1, and A2 information for each chromosome
  fwrite(ref[ref$CHR == chr_i, c('SNP','A1','A2'), with=F], paste0(tmp_dir,'/map.txt'), row.names=F, quote=F, sep=' ')

  # Paste files together in batches
  # Set number of batches according to the number of score files to combine
  num_batches <- max(c(1, min(c(opt$n_cores, floor(nrow(score_files) / 2)))))
  tmp_score_files <- paste0(tmp_dir, '/score_shard.', seq_len(nrow(score_files)), '.', chr_i, '.txt')
  set.seed(1)
  batches <- split(sample(tmp_score_files), rep(1:num_batches, length.out = length(tmp_score_files)))
  log_add(log_file = log_file, message = paste0('Aggregating score files in ', num_batches,' batches.'))
  foreach(i = 1:length(batches), .combine = c, .options.multicore = list(preschedule = FALSE)) %dopar% {
    join_score_files(batches[[i]], paste0(tmp_dir,'/tmp_batch_',i))
    unlink(batches[[i]])
  }

  # Paste batches together
  log_add(log_file = log_file, message = paste0('Aggregating batched score files.'))
  tmp_batch_files <- paste0(tmp_dir,'/tmp_batch_',1:length(batches))
  join_score_files(c(paste0(tmp_dir,'/map.txt'), tmp_batch_files), paste0(tmp_dir, '/all_score.txt'))
  unlink(tmp_batch_files)

  # Perform polygenic risk scoring
  scores_i <-
    plink_score(
      pfile = opt$ref_plink_chr,
      chr = chr_i,
      plink2 = opt$plink2,
      score = paste0(tmp_dir,'/all_score.txt'),
      threads = opt$n_cores
    )
  
  # Sum scores across chromosomes
  if(chr_i == CHROMS[1]){
    scores_ids <- scores_i[, 1:2, with = F]
    current_scores <- as.matrix(scores_i[, -1:-2, with = FALSE])
    scores <- current_scores
  } else {
    check_sample_order(scores_ids, scores_i)
    if(!identical(colnames(scores), names(scores_i)[-c(1, 2)])) stop('Chromosome score columns differ.')
    current_scores <- as.matrix(scores_i[, -1:-2, with = FALSE])
    scores <- scores + current_scores
  }

  unlink(paste0(tmp_dir, c('/all_score.txt', '/map.txt')))
  rm(scores_i)
  rm(current_scores)
  gc()
}

# Combine score with IDs
scores<-data.table(scores_ids,
                   scores)

###
# Scale the polygenic scores based on the reference
###

log_add(log_file = log_file, message = paste0('Adjusting PGS for ancestry.'))

pop_data <- read_pop_data(paste0(refdir, '/ref.pop.txt'))
validate_sample_ids(pop_data, 'reference population data')

for(i in 1:nrow(score_files)){
  scores_i <- scores[, c('FID','IID', names(scores)[grepl(paste0('^score_file_', i, '\\.'), names(scores))]), with=F]
  names(scores_i) <- gsub(paste0('^score_file_', i, '\\.'), 'SCORE_', names(scores_i))

  output_i <- paste0(outdir, '/reference/pgs_score_files/', score_files$method[i], '/', score_files$name[i], '/ref-', score_files$name[i])

  if(opt$continuous){
    # Derive trans-ancestry PGS models and estimate PGS residual scale
    model_trans_pgs(scores=scores_i, pcs=paste0(outdir, '/reference/pc_score_files/TRANS/ref-TRANS-pcs.profiles'), output=output_i)
  }
  
  # Calculate scale within each reference population
  for(pop_i in unique(pop_data$POP)){
    unlink(paste0(output_i, '-', pop_i, c('.scale', '.profiles')))
    ref_pgs_scale_i <- score_mean_sd(scores = scores_i, keep = pop_data[pop_data$POP == pop_i, c('FID','IID'), with=F])
    fwrite(ref_pgs_scale_i, paste0(output_i, '-', pop_i, '.scale'), row.names = F, quote=F, sep=' ', na='NA')
    scores_i_pop <- scores_i[pop_data[POP == pop_i, .(FID, IID)], on = .(FID, IID), nomatch = 0]
    scores_i_pop<-score_scale(score=scores_i_pop, ref_scale=ref_pgs_scale_i)
    fwrite(scores_i_pop, paste0(output_i, '-', pop_i, '.profiles'), sep=' ', na='NA', quote=F)
  }
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
