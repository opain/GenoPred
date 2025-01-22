#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NA, type='character',
              help="Path to per chromosome reference PLINK files [required]"),
  make_option("--ref_freq_chr", action="store", default=NULL, type='character',
              help="Path to per chromosome reference PLINK2 .afreq files [required]"),
  make_option("--ref_pcs", action="store", default=NULL, type='character',
              help="Reference PCs for continuous ancestry correction [optional]"),
  make_option("--pop_data", action="store", default=NULL, type='character',
              help="File containing the population code and location of the keep file [required]"),
  make_option("--plink", action="store", default='plink', type='character',
              help="Path PLINK v1.9 software binary [optional]"),
  make_option("--plink2", action="store", default='plink2', type='character',
              help="Path PLINK v2 software binary [optional]"),
  make_option("--output", action="store", default='NA', type='character',
              help="Path for output files [required]"),
  make_option("--memory", action="store", default=5000, type='numeric',
              help="Memory limit [optional]"),
  make_option("--sumstats", action="store", default=NULL, type='character',
              help="Comma-seperated list of GWAS summary statistics [required]"),
  make_option("--populations", action="store", default=NULL, type='character',
              help="Comma-seperated list of population codes matching GWAS [required]"),
  make_option("--ldak", action="store", default=NA, type='character',
              help="Path to ldak v5.2 executable [required]"),
  make_option("--quickprs_ldref", action="store", default=NA, type='character',
              help="Path to folder containing ldak quickprs reference [required]"),
  make_option("--quickprs_multi_ldref", action="store", default=NA, type='character',
              help="Path to folder containing ldak quickprs_multi reference [required]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
              help="Number of cores for parallel computing [optional]"),
  make_option("--prs_model", action="store", default='bayesr', type='character',
              help="Model used for deriving SNP-weights [optional]"),
  make_option("--genomic_control", action="store", default=F, type='logical',
              help="Logical indicating whether genomic control was applied to GWAS [optional]"),
  make_option("--xwing_repo", action="store", default=NULL, type='character',
              help="Path to X-WING repo [required]"),
  make_option("--test", action="store", default=NA, type='character',
              help="Specify number of SNPs to include [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')

# Check required inputs
if(is.null(opt$ref_plink_chr)){
  stop('--ref_plink_chr must be specified.\n')
}
if(is.null(opt$ref_freq_chr)){
  stop('--ref_freq_chr must be specified.\n')
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
if(is.null(opt$ldak)){
  stop('--ldak must be specified.\n')
}
if(is.null(opt$quickprs_ldref)){
  stop('--quickprs_ldref must be specified.\n')
}
if(is.null(opt$quickprs_multi_ldref)){
  stop('--quickprs_multi_ldref must be specified.\n')
}
if(is.null(opt$populations)){
  stop('--populations must be specified.\n')
}

# Create output directory
opt$output_dir <- paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'quickprs_multi.R', start.time = start.time)

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

# Split opt$populations
populations<-unlist(strsplit(opt$populations, ','))

#####
# Format the sumstats
#####

gwas_N<-NULL
for(i in 1:length(sumstats)){
  log_add(log_file = log_file, message = paste0('Reading in GWAS ', i))

  # Read in, check and format GWAS summary statistics
  gwas <- read_sumstats(sumstats = sumstats[i], chr = CHROMS, log_file = log_file, req_cols = c('CHR','BP','SNP','A1','A2','BETA','SE','N','P'))

  # Format for LDAK
  snplist <- gwas$SNP
  gwas$Z <- gwas$BETA / gwas$SE
  gwas$Predictor<-paste0(gwas$CHR, ':', gwas$BP)
  gwas<-gwas[,c('Predictor','A1','A2','N','Z','CHR','BP','Predictor','BETA','P')]
  names(gwas)<-c('Predictor','A1','A2','n','Z','CHR','BP','SNP','BETA','P')
  gwas_N <- c(gwas_N, round(mean(gwas$n), 0))

  # Check overlap between GWAS and LDAK reference
  quickprs_ldref_pop_i <- paste0(opt$quickprs_ldref, '/', populations[i])
  ldak_hm3_file <- list.files(quickprs_ldref_pop_i)
  ldak_hm3_file <- ldak_hm3_file[grepl('.cors.bim', ldak_hm3_file)][1]
  ldak_hm3 <- fread(paste0(quickprs_ldref_pop_i, '/', ldak_hm3_file))
  ldak_hm3 <- ldak_hm3[ldak_hm3$V1 %in% CHROMS,]
  ref_overlap <- sum(gwas$Predictor %in% ldak_hm3$V2) / nrow(ldak_hm3)

  log_add(log_file = log_file, message = paste0('GWAS ', i,'-reference overlap is ', round(ref_overlap * 100, 2), '%.'))

  # Subset GWAS to LDAK reference data
  gwas <- gwas[gwas$Predictor %in% ldak_hm3$V2, ]

  # Output formatted sumstats
  fwrite(gwas, paste0(tmp_dir,'/GWAS_sumstats_temp', i, '.txt'), sep=' ')
}

# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

#####
# Run QuickPRS using full sumstats
#####

log_add(log_file = log_file, message = 'Running QuickPRS using full sumstats.')

score_full <- list()
for(i in 1:length(sumstats)){
  score_full[[populations[i]]] <- quick_prs(
    sumstats = paste0(tmp_dir,'/GWAS_sumstats_temp', i, '.txt'),
    ref_dir = paste0(opt$quickprs_ldref, '/', populations[i]),
    genomic_control = opt$genomic_control,
    n_cores = opt$n_cores,
    prs_model = opt$prs_model)
}

#####
# Subsample GWAS sumstats
#####

log_add(log_file = log_file, message = 'Subsampling sumstats.')

dir.create(paste0(tmp_dir,'/LEOPARD/sampled_sumstats'), recursive = T)

for(i in 1:length(sumstats)){
  quickprs_multi_ldref_pop_i <- paste0(opt$quickprs_multi_ldref, '/', populations[i])
  ref_files <- list.files(quickprs_multi_ldref_pop_i)
  ref_files <- gsub('.bed', '', ref_files[grepl('subset_1.bed', ref_files)])
  system(
    paste0(
      'Rscript ',
      opt$xwing_repo,
      '/LEOPARD_Sim.R ',
      '--sumstats ', tmp_dir, '/GWAS_sumstats_temp', i, '.txt ',
      '--n_gwas ', gwas_N[i], ' ',
      '--train_prop 0.75 ',
      '--ref_prefix ', quickprs_multi_ldref_pop_i, '/', ref_files, ' ',
      '--seed 1 ',
      '--rep 4 ',
      '--out_prefix ', tmp_dir, '/LEOPARD/sampled_sumstats/GWAS_', i
    )
  )
}

#####
# Run QuickPRS using subsampled sumstats
#####

log_add(log_file = log_file, message = 'Running QuickPRS on subsampled sumstats.')

score_subset <- list()
for(i in 1:length(sumstats)){
  dir.create(paste0(tmp_dir, '/LEOPARD/post_targ_', populations[i]))
  score_subset[[populations[i]]] <- list()
  for(j in 1:4){
    sumstats_tmp <- fread(paste0(tmp_dir, '/LEOPARD/sampled_sumstats/GWAS_', i, '_rep', j, '_train.txt'))
    n_tmp <- fread(paste0(tmp_dir, '/LEOPARD/sampled_sumstats/GWAS_', i, '_rep', j, '_train_valid_N.txt'))
    sumstats_tmp$n <- n_tmp$N_train
    sumstats_tmp$Direction <-sign(sumstats_tmp$BETA)
    sumstats_tmp <- sumstats_tmp[, c('Predictor','A1','A2','n','P','Direction'), with = F]
    fwrite(sumstats_tmp, paste0(tmp_dir, '/LEOPARD/sampled_sumstats/GWAS_', i, '_rep', j, '_train.reformat.txt'), sep=' ')

    score_subset[[populations[i]]][[paste0('subset_', j)]] <- quickprs(
      sumstats = paste0(tmp_dir, '/LEOPARD/sampled_sumstats/GWAS_', i, '_rep', j, '_train.reformat.txt'),
      quickprs_ldref = paste0(opt$quickprs_ldref, '/', populations[i]),
      quickprs_multi_ldref = paste0(opt$quickprs_multi_ldref, '/', populations[i]),
      genomic_control = opt$genomic_control,
      n_cores = opt$n_cores,
      ref_subset = '2',
      prs_model = opt$prs_model)

    fwrite(score_subset[[populations[i]]][[paste0('subset_', j)]],
           paste0(tmp_dir, '/LEOPARD/post_targ_', populations[i], '/output_', j, '_', populations[i], '_Post.txt'),
           sep=' ')
  }
}

#####
# Estimating the linear combination weights
#####

log_add(log_file = log_file, message = 'Estimating the linear combination weights.')

for(targ_pop in populations){
  dir.create(paste0(tmp_dir,'/LEOPARD/weights_', targ_pop), recursive = T)
  quickprs_multi_ldref_pop_i <- paste0(opt$quickprs_multi_ldref, '/', targ_pop)
  ref_files <- list.files(quickprs_multi_ldref_pop_i)
  ref_files <- gsub('.bed', '', ref_files[grepl('subset_3.bed', ref_files)])

  for(j in 1:4){
    targ_gwas_valid_n<-fread(paste0(tmp_dir, '/LEOPARD/sampled_sumstats/GWAS_', which(populations == targ_pop), '_rep', j, '_train_valid_N.txt'))$N_valid

    sumstats_tmp <- fread(paste0(tmp_dir, '/LEOPARD/sampled_sumstats/GWAS_', which(populations == targ_pop), '_rep', j, '_valid.txt'))
    sumstats_tmp$SNP<-sumstats_tmp$Predictor
    sumstats_tmp <- sumstats_tmp[, c('SNP','CHR','BP','A1','A2','BETA','P'), with = F]
    fwrite(sumstats_tmp, paste0(tmp_dir, '/LEOPARD/sampled_sumstats/GWAS_', which(populations == targ_pop), '_rep', j, '_valid.reformat.txt'), sep=' ')

    system(paste0(
      'Rscript ', opt$xwing_repo, '/LEOPARD_Weights.R ',
      '--beta_file ', tmp_dir, '/LEOPARD/post_targ_', populations[1], '/output_', j, '_', populations[1], '_Post.txt',',',tmp_dir, '/LEOPARD/post_targ_', populations[2], '/output_', j, '_', populations[2], '_Post.txt ',
      '--valid_file ', tmp_dir, '/LEOPARD/sampled_sumstats/GWAS_', which(populations == targ_pop), '_rep', j, '_valid.reformat.txt ',
      '--n_valid ', targ_gwas_valid_n ,' ',
      '--ref_prefix ', quickprs_multi_ldref_pop_i, '/', ref_files, ' ',
      '--out ', tmp_dir,'/LEOPARD/weights_', targ_pop,'/output_LEOPARD_weights_rep', j, '.txt'
    ))
  }
}

# Average weights across repeats
mix_weights <- calculate_avg_weights(populations = populations, leopard_dir = paste0(tmp_dir,'/LEOPARD'), log_file = log_file)

####
# Combine score files
####

log_add(log_file = log_file, message = 'Creating score file.')

# Combine the scores from each population
score_all <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c('SNP','A1','A2'), all = TRUE), score_full)
names(score_all)[grepl('BETA', names(score_all))]<-paste0('SCORE_targ_', populations)
score_all[is.na(score_all)]<-0

# Read in reference SNP data
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)

# Change IDs to RSIDs
ref$Predictor<-paste0(ref$CHR,':',ref$BP)
names(score_all)[names(score_all) == 'SNP'] <- 'Predictor'
score_all<-merge(score_all, ref[,c('Predictor','SNP'), with=F], by='Predictor')
score_all<-score_all[, c('SNP', 'A1', 'A2', names(score_all)[grepl('SCORE_', names(score_all))]), with=F]

# Flip effects to match reference alleles
score_all <- map_score(ref = ref, score = score_all)

# Calculate linear combination of scores using mixing weights for each target population
score_weighted <- score_all
for(targ_pop in populations){
  # Read in the .freq file for target population
  freq_data <- read_frq(freq_dir = opt$ref_freq_chr, population = targ_pop, chr = CHROMS)
  
  # Centre SNP-weights for target population
  score_i <- centre_weights(score = score_all, freq = freq_data, ref = ref)
  
  # Linearly combine scores using mixing weights for target population
  score_weighted[[paste0('SCORE_targ_', targ_pop, '_weighted')]] <-
    calculate_weighted_scores(score = score_i,
                              mix_weights = mix_weights,
                              targ_pop = targ_pop)
}

# Reduce number of significant figures to save space
score_weighted[, (4:ncol(score_all)) := lapply(.SD, signif, digits = 7), .SDcols = 4:ncol(score_all)]

fwrite(score_weighted, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

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

if(!is.null(opt$ref_pcs)){
  log_add(log_file = log_file, message = 'Deriving trans-ancestry PGS models...')
  # Derive trans-ancestry PGS models and estimate PGS residual scale
  model_trans_pgs(scores=ref_pgs, pcs=opt$ref_pcs, output=opt$output)
}

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
