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

# Check whether score files or target genetic data are newer than target pgs
if(!is.null(score_files)){
  ref_pcs_file_time <- NULL
  if(opt$continuous){
    ref_pcs_file<-paste0(outdir, '/reference/pc_score_files/TRANS/ref-TRANS-pcs.EUR.scale')
    ref_pcs_file_time <- file.info(ref_pcs_file)$mtime
  }
  
  score_files_to_do <- data.table()
  for(i in 1:nrow(score_files)){
    pgs_i <- paste0(outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i], '-EUR.profiles')
    score_i <- paste0(outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i], '.score.gz')
    if(!file.exists(pgs_i)){
      score_files_to_do <- rbind(score_files_to_do, score_files[i,])
    } else {
      score_i_time <- file.info(score_i)$mtime
      pgs_i_time <- file.info(pgs_i)$mtime
      if (score_i_time > pgs_i_time | (!is.null(ref_pcs_file_time) && ref_pcs_file_time > pgs_i_time)) {
        score_files_to_do <- rbind(score_files_to_do, score_files[i,])
        system(paste0('rm ', pgs_i))
      }
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

# Read in reference SNP data
refdir <- read_param(config = opt$config, param = 'refdir', return_obj = F)
if(as.logical(read_param(config = opt$config, param = 'restrict_to_target_variants', return_obj = F)) |
   as.logical(read_param(config = opt$config, param = 'dense_reference', return_obj = F))){
  opt$ref_plink_chr <- paste0(outdir, '/reference/ref/ref.chr')
} else {
  opt$ref_plink_chr <- paste0(refdir, '/ref.chr')
}

# Read in reference SNP data
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('CHR','SNP','A1','A2'), with=F]

# We will process score files and perform target scoring for one chromosome for efficiency
score_files_orig <- score_files
for(chr_i in CHROMS){
  score_files <- score_files_orig
  log_add(log_file = log_file, message = '########################')
  log_add(log_file = log_file, message = paste0('Processing chromosome ', chr_i,':'))

  # Initialize a matrix running PGS sum
  if(chr_i == CHROMS[1]){
    
    scores_ids <- fread(cmd=paste0('cut -f 1 ', opt$ref_plink_chr, CHROMS[1], '.psam'), header = T)
    names(scores_ids)<-gsub('\\#', '', names(scores_ids))
    scores_ids <- data.table(FID = scores_ids$IID, IID = scores_ids$IID)
    
    cols <- NULL
    for(i in 1:nrow(score_files)){
      cols_i <- names(fread(cmd = paste0('zcat ', outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i],'.score.gz | head -n 1'), header=T))[-1:-3]
      cols_i <- gsub('SCORE_', paste0('score_file_', i,'.'), cols_i)
      cols <- c(cols, cols_i)
    }
    
    scores <- matrix(
      nrow = nrow(scores_ids),
      ncol = length(cols),
      data = 0
    )
  }
  
  #####
  # Combine score files
  #####
  # Identify score files with data for chromosome
  check <- foreach(i = 1:nrow(score_files), .combine = rbind, .options.multicore = list(preschedule = FALSE)) %dopar% {
    if(score_files$method[i] == 'external'){
      ss_file <- paste0(outdir,'/reference/pgs_score_files/external/',score_files$name[i],'/ref-',score_files$name[i],'.unmapped.score.gz')
    } 
    if(score_files$method[i] != 'external') {
      if(!(score_files$method[i] %in% c('prscsx','xwing')) & !grepl('_multi|tlprs_', score_files$method[i])){
        ss_file <- paste0(outdir,'/reference/gwas_sumstat/',score_files$name[i],'/',score_files$name[i],'-cleaned.gz')
      } else {
        gwas_groups <- read_param(config = opt$config, param = 'gwas_groups', return_obj = T)
        gwas_groups <- gwas_groups[gwas_groups$name == score_files$name[i],]
        gwas_groups <- unlist(strsplit(gwas_groups$gwas, ','))
        ss_file <- paste0(outdir,'/reference/gwas_sumstat/',gwas_groups[1],'/',gwas_groups[1],'-cleaned.gz')
      }
    }
    ss_names <- fread(cmd = paste0('zcat ', ss_file,' | head -n 1'), nrows = 0, header = T)
    exit_code <- system(paste0('zcat ', ss_file, " | awk 'NR>1 && $", which(names(ss_names) == 'CHR'), "+0==", chr_i, " { found=1; exit } END {if (!found) exit 1}'"))
    
    data.frame(i = i, zero_only = exit_code == 1)
  }
  score_files <- score_files[check$i[!check$zero_only],]

  if(nrow(score_files) > 0){
    # Create row number index to subset score files by chromosome
    row_index <- format(which(ref$CHR == chr_i) + 1, scientific = FALSE)
    write.table(row_index, paste0(tmp_dir,'/row_index.txt'), row.names=F, quote=F, col.names = F)
  
    # Create file containing SNP, A1, and A2 information for each chromosome
    fwrite(ref[ref$CHR == chr_i, c('SNP','A1','A2'), with=F], paste0(tmp_dir,'/map.txt'), row.names=F, quote=F, sep=' ')
  
    # Extract process score files for each name (gwas/score) in parallel
    foreach(i = 1:nrow(score_files), .combine = c, .options.multicore = list(preschedule = FALSE)) %dopar% {
      system(paste0(
        'zcat ', outdir, '/reference/pgs_score_files/', score_files$method[i],'/', score_files$name[i],'/ref-',score_files$name[i],".score.gz | ",
        'awk \'NR==FNR {rows[$1]; next} FNR==1 || FNR in rows\' ', paste0(tmp_dir,'/row_index.txt'), ' - | ',  # Corrected to retain the header and process indexed rows
        "cut -d' ' --complement -f1-3 | ",  # Keep relevant columns, remove first 3
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
        pfile = opt$ref_plink_chr,
        chr = chr_i,
        plink2 = opt$plink2,
        score = paste0(tmp_dir,'/all_score.txt'),
        threads = opt$n_cores
      )
  
    # Sum scores across chromosomes
    # Reorder scores to match matrix
    match_idx <- match(paste(scores_ids$FID, scores_ids$IID),
                       paste(scores_i$FID, scores_i$IID))
    
    # In-place addition: for each score column
    scores_i<-as.matrix(scores_i[,-1:-2])
    for (j in cols) {
      scores[, which(cols == j)] <- scores[, which(cols == j)] + scores_i[match_idx, which(colnames(scores_i) == j)]
    }

    system(paste0('rm ', tmp_dir, '/all_score.txt'))
    system(paste0('rm ', tmp_dir, '/row_index.txt'))
    system(paste0('rm ', tmp_dir, '/map.txt'))
    rm(scores_i)
    gc()
  }
}

# Combine score with IDs
scores <- data.table(scores)
setnames(scores, cols)
scores <- cbind(scores_ids, scores)

###
# Scale the polygenic scores based on the reference
###

log_add(log_file = log_file, message = paste0('Adjusting PGS for ancestry.'))

pop_data <- read_pop_data(paste0(refdir, '/ref.pop.txt'))

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
    ref_pgs_scale_i <- score_mean_sd(scores = scores_i, keep = pop_data[pop_data$POP == pop_i, c('FID','IID'), with=F])
    fwrite(ref_pgs_scale_i, paste0(output_i, '-', pop_i, '.scale'), row.names = F, quote=F, sep=' ', na='NA')
    scores_i_pop<-scores_i[paste0(scores_i$FID, '_', scores_i$IID) %in% paste0(pop_data$FID[pop_data$POP == pop_i], '_', pop_data$IID[pop_data$POP == pop_i]),]
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
