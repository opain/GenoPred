#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NULL, type='character',
      help="Path to per chromosome reference PLINK files [required]"),
  make_option("--xwing_repo", action="store", default=NULL, type='character',
      help="Path to X-WING repo [required]"),
  make_option("--logodetect_ref", action="store", default=NULL, type='character',
      help="Path LOGODetect reference data [required]"),
  make_option("--panther_ref", action="store", default=NULL, type='character',
      help="Path PANTHER reference data [required]"),
  make_option("--leopard_ref", action="store", default=NULL, type='character',
      help="Path LEOPARD reference data [required]"),
  make_option("--panther_leopard_ref", action="store", default=NULL, type='character',
      help="Path subsampled reference data for PANTHER/LEOPARD [required]"),
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
  make_option("--n_cores", action="store", default=1, type='numeric',
      help="Number of cores for parallel computing [optional]"),
  make_option("--test", action="store", default=NA, type='character',
      help="Specify number of SNPs to include [optional]"),
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
if(is.null(opt$xwing_repo)){
  stop('--xwing_repo must be specified.\n')
}
if(is.null(opt$logodetect_ref)){
  stop('--logodetect_ref must be specified.\n')
}
if(is.null(opt$panther_ref)){
  stop('--panther_ref must be specified.\n')
}
if(is.null(opt$leopard_ref)){
  stop('--leopard_ref must be specified.\n')
}
if(is.null(opt$panther_leopard_ref)){
  stop('--panther_leopard_ref must be specified.\n')
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

# Split opt$sumstats
populations<-unlist(strsplit(opt$populations, ','))

gwas_N<-NULL
for(i in 1:length(sumstats)){

  #####
  # Read in sumstats
  #####

  log_add(log_file = log_file, message = 'Reading in GWAS.')

  # Read in, check and format GWAS summary statistics
  gwas <- read_sumstats(sumstats = sumstats[i], chr = CHROMS, log_file = log_file, req_cols = c('CHR', 'SNP', 'BP', 'A1', 'A2', 'BETA', 'P', 'N'))

  # Store average sample size
  gwas_N <- c(gwas_N, round(mean(gwas$N), 0))
  gwas$N<-NULL

  fwrite(gwas, paste0(tmp_dir, '/GWAS_sumstats_',i,'_temp.txt'), sep=' ')

  rm(gwas)
  gc()


}

# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

#####
# Process sumstats using X-WING
#####

# Make a list of GWAS pairs to be analysed
gwas_pairs <- combn(1:length(sumstats), 2)
gwas_pairs <- as.matrix(data.frame(gwas_1 = gwas_pairs[1, ], gwas_2 = gwas_pairs[2, ]))

# Run LOGODetect
dir.create(paste0(tmp_dir,'/LOGODetect'))

for(i in 1:nrow(gwas_pairs)){
  system(paste0(
    'Rscript ', opt$xwing_repo, '/LOGODetect.R --sumstats ', paste0(paste0(tmp_dir, '/GWAS_sumstats_', gwas_pairs[i,],'_temp.txt'), collapse=','),' --n_gwas ', paste0(gwas_N[gwas_pairs[i,]], collapse=','), ' --ref_dir ', opt$logodetect_ref,' --pop ', paste0(populations[gwas_pairs[i,]], collapse=',') ,' --block_partition ', opt$xwing_repo,'/block_partition.txt --gc_snp ', opt$xwing_repo,'/1kg_hm3_snp.txt --out_dir ', tmp_dir, '/LOGODetect --n_cores ', opt$n_core,' --target_pop ', opt$populations,' --n_topregion 1000'
  ))
}

##
# Run PANTHER
##
# Create a temporary reference bim files for X-WING to match
pvar <- read_pvar(opt$ref_plink_chr, chr = CHROMS)
pvar$POS<-0
for(i in CHROMS){
  write.table(pvar[pvar$CHR == i, c('CHR','SNP','POS','BP','A1','A2'), with=F], paste0(tmp_dir,'/ref.chr',i,'.bim'), col.names=F, row.names=F, quote=F)
}

rm(pvar)
gc()

combinations <- expand.grid(targ_pop = populations, pst_pop = populations, chr = CHROMS)

file.remove(paste0(tmp_dir, '/checker.txt'))
log <- foreach(i = 1:nrow(combinations), .combine = c, .options.multicore = list(preschedule = FALSE)) %dopar% {
  if(!file.exists(paste0(tmp_dir, '/checker.txt'))) {
    targ_pop <- combinations$targ_pop[i]
    pst_pop <- combinations$pst_pop[i]
    chr <- combinations$chr[i]

    # Create directories
    dir.create(paste0(tmp_dir, '/PANTHER/post_targ_', targ_pop), recursive = TRUE)
    dir.create(paste0(tmp_dir, '/PANTHER/post_collect_targ_', targ_pop), recursive = TRUE)

    command <- paste0(
      'python ', opt$xwing_repo, '/PANTHER.py ',
      '--ref_dir ', opt$panther_ref, ' ',
      '--bim_prefix ', tmp_dir,'/ref.chr',chr,' ',
      '--sumstats ', paste0(paste0(tmp_dir, '/GWAS_sumstats_', 1:length(sumstats),'_temp.txt'), collapse=','), ' ',
      '--n_gwas ', paste(gwas_N, collapse=','), ' ',
      '--anno_file ', paste0(paste0(tmp_dir, '/LOGODetect/targ_', targ_pop, '_annot_', populations, '.txt'), collapse=','), ' ',
      '--chrom ', chr, ' ',
      '--pop ', opt$populations ,' ',
      '--target_pop ', targ_pop,' ',
      '--pst_pop ', pst_pop, ' ',
      '--out_name output ',
      '--seed 1 ',
      '--out_dir ', tmp_dir, '/PANTHER/post_targ_', targ_pop
    )

    # Run command
    log_i <- system(command)

    # Check for an error
    if(log_i != 0){
      write("", paste0(tmp_dir, '/checker.txt'))
    }
  }
}

##
# Run LEOPARD to subsample GWAS
##
dir.create(paste0(tmp_dir,'/LEOPARD/sampled_sumstats'), recursive = T)

for(i in 1:length(sumstats)){
  system(paste0(
    'Rscript ', opt$xwing_repo, '/LEOPARD_Sim.R ',
    '--sumstats ', tmp_dir, '/GWAS_sumstats_', i,'_temp.txt ',
    '--n_gwas ', gwas_N[i], ' ',
    '--train_prop 0.75 ',
    '--ref_prefix ', opt$leopard_ref,'/', populations[i], '/', populations[i], '_part1 ',
    '--seed ', opt$seed, ' ',
    '--rep 4 ',
    '--out_prefix ', tmp_dir,'/LEOPARD/sampled_sumstats/GWAS_', i
  ))
}

##
# Run PANTHER on subsampled GWAS
##

combinations <- expand.grid(targ_pop = populations, pst_pop = populations, chr = CHROMS, index = 1:4)

file.remove(paste0(tmp_dir, '/checker.txt'))
log <- foreach(i = 1:nrow(combinations), .combine = c, .options.multicore = list(preschedule = FALSE)) %dopar% {
  if(!file.exists(paste0(tmp_dir, '/checker.txt'))) {
    targ_pop <- combinations$targ_pop[i]
    pst_pop <- combinations$pst_pop[i]
    chr <- combinations$chr[i]
    index <- combinations$index[i]

    dir.create(paste0(tmp_dir,'/LEOPARD/post_targ_', targ_pop), recursive = T)
    dir.create(paste0(tmp_dir,'/LEOPARD/post_collect_targ_', targ_pop), recursive = T)

    sumstats_i <- paste0(tmp_dir, '/GWAS_sumstats_', 1:length(sumstats),'_temp.txt')
    sumstats_i[populations == targ_pop] <- paste0(tmp_dir,'/LEOPARD/sampled_sumstats/GWAS_', which(populations == targ_pop), '_rep', index, '_train.txt')

    targ_gwas_train_n<-fread(paste0(tmp_dir,'/LEOPARD/sampled_sumstats/GWAS_', which(populations == targ_pop), '_rep', index, '_train_valid_N.txt'))$N_train

    gwas_N_i<-gwas_N
    gwas_N_i[populations == targ_pop] <- targ_gwas_train_n

    command<-paste0(
      'python ', opt$xwing_repo, '/PANTHER.py ',
      '--ref_dir ', opt$panther_leopard_ref, ' ',
      '--bim_prefix ', tmp_dir,'/ref.chr',chr,' ',
      '--sumstats ', paste0(sumstats_i, collapse=','), ' ',
      '--n_gwas ', paste(gwas_N_i, collapse=','), ' ',
      '--anno_file ', paste0(paste0(tmp_dir, '/LOGODetect/targ_', targ_pop, '_annot_', populations, '.txt'), collapse=','), ' ',
      '--chrom ', chr, ' ',
      '--pop ', opt$populations ,' ',
      '--target_pop ', targ_pop,' ',
      '--pst_pop ', pst_pop, ' ',
      '--out_name output_', index, ' ',
      '--seed 1 ',
      '--out_dir ', tmp_dir, '/LEOPARD/post_targ_', targ_pop
    )

    # Run command
    log_i <- system(command)

    # Check for an error
    if(log_i != 0){
      write("", paste0(tmp_dir, '/checker.txt'))
    }
  }
}

for(targ_pop in populations){
  for(pst_pop in populations){
    for(i in 1:4){
      system(paste0('cat ', tmp_dir, '/LEOPARD/post_targ_', targ_pop, '/output_', i,'_', pst_pop, '_pst_eff_phiauto_chr*.txt > ', tmp_dir, '/LEOPARD/post_targ_', targ_pop, '/output_', i,'_', pst_pop, '_Post.txt'))
      system(paste0("sed  -i '1iCHR\tSNP\tBP\tA1\tA2\tBETA' ", tmp_dir, '/LEOPARD/post_targ_', targ_pop, '/output_', i, '_', pst_pop, '_Post.txt'))
    }
  }
}

##
# Run LEOPARD to to find optimal weights for GWAS from each population
##

# Estimating the linear combination weights
for(targ_pop in populations){
  dir.create(paste0(tmp_dir,'/LEOPARD/weights_', targ_pop), recursive = T)

  for(i in 1:4){

    targ_gwas_valid_n<-fread(paste0(tmp_dir,'/LEOPARD/sampled_sumstats/GWAS_', which(populations == targ_pop), '_rep', i, '_train_valid_N.txt'))$N_valid

    system(paste0(
      'Rscript ', opt$xwing_repo, '/LEOPARD_Weights.R --beta_file ', paste0(paste0(tmp_dir, '/LEOPARD/post_targ_', targ_pop, '/output_', i, '_', populations, '_Post.txt'), collapse=','), ' --valid_file ', tmp_dir,'/LEOPARD/sampled_sumstats/GWAS_', which(populations == targ_pop), '_rep', i, '_valid.txt --n_valid ', targ_gwas_valid_n ,' --ref_prefix ', opt$leopard_ref,'/', targ_pop, '/', targ_pop, '_part3 --out ', tmp_dir,'/LEOPARD/weights_', targ_pop,'/output_LEOPARD_weights_rep', i, '.txt'
    ))
  }
}

# Average weights across repeats
for(targ_pop in populations){
  system(paste0(
    'Rscript ', opt$xwing_repo, '/LEOPARD_Avg.R --weights_prefix ', tmp_dir,'/LEOPARD/weights_', targ_pop,'/output_LEOPARD_weights_rep --rep 4 --out ', tmp_dir,'/LEOPARD/weights_', targ_pop,'/output_LEOPARD_weights_avg.txt'
  ))
}

####
# Combine score files
####

# We should combine the raw PANTHER score files for each population,
# and then combine using mixing weights for each population
score_all<-NULL
for(targ_pop in populations){
  mix_weights<-fread(paste0(tmp_dir,'/LEOPARD/weights_', targ_pop,'/output_LEOPARD_weights_avg.txt'))
  score_pop <- NULL
  for(pst_pop in populations){
    score_i<-NULL
    for(chr in CHROMS){
      score_i_chr<-fread(paste0(tmp_dir, '/PANTHER/post_targ_', targ_pop, '/output_', pst_pop, '_pst_eff_phiauto_chr', chr,'.txt'), header=F)
      score_i<-rbind(score_i, score_i_chr)
    }

    score_i$V6 <- score_i$V6 * mix_weights$Weights[grepl(paste0(pst_pop,'_Post.txt'), mix_weights$Path)]

    names(score_i)<-c('CHR','SNP', 'BP', 'A1', 'A2', paste0('SCORE_targ_', targ_pop, '_pst_', pst_pop))
    score_i<-score_i[, c('SNP', 'A1', 'A2', paste0('SCORE_targ_', targ_pop, '_pst_', pst_pop)), with=F]

    if(is.null(score_pop)){
      score_pop<-score_i
    } else {
      score_pop<-merge(score_pop, score_i, by=c('SNP','A1','A2'), all=T)
    }
  }
  # Take average of weighted scores
  score_pop[is.na(score_pop)]<-0
  score_pop[[paste0('SCORE_targ_', targ_pop, '_weighted')]] <- rowSums(score_pop[, grepl('SCORE_', names(score_pop)), with = F])

  if(is.null(score_all)){
    score_all<-score_pop
  } else {
    score_all<-merge(score_all, score_pop, by=c('SNP','A1','A2'), all=T)
  }
}

# Flip effects to match reference alleles
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]
score_new <- map_score(ref = ref, score = score_all)

# Reduce number of significant figures to save space
score_new[, (4:ncol(score_new)) := lapply(.SD, signif, digits = 7), .SDcols = 4:ncol(score_new)]

fwrite(score_new, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

if(file.exists(paste0(opt$output,'.score.gz'))){
  system(paste0('rm ',opt$output,'.score.gz'))
}

system(paste0('gzip ',opt$output,'.score'))

# Record end time of test
if(!is.na(opt$test)){
  test_finish(log_file = log_file, test_start.time = test_start.time)
}

####
# Calculate mean and sd of polygenic scores
####

log_add(log_file = log_file, message = 'Calculating polygenic scores in reference.')

# Calculate scores in the full reference
ref_pgs <- plink_score(pfile = opt$ref_plink_chr, chr = CHROMS, plink2 = opt$plink2, score = paste0(opt$output,'.score.gz'), threads = opt$n_cores)

# Calculate scale within each reference population
pop_data <- read_pop_data(opt$pop_data)

for(pop_i in unique(pop_data$POP)){
  ref_pgs_scale_i <- score_mean_sd(scores = ref_pgs, keep = pop_data[pop_data$POP == pop_i, c('FID','IID'), with=F])
  fwrite(ref_pgs_scale_i, paste0(opt$output, '-', pop_i, '.scale'), row.names = F, quote=F, sep=' ', na='NA')
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()
