#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--config", action="store", default=NA, type='character',
              help="Path to config file [required]"),
  make_option("--method", action="store", default=NA, type='character',
              help="PGS method [required]"),
  make_option("--gwas_group", action="store", default=NA, type='character',
              help="GWAS group [required]"),
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
if(is.null(opt$config)){
  stop('--config must be specified.\n')
}
if(is.null(opt$method)){
  stop('--method must be specified.\n')
}
if(is.null(opt$gwas_group)){
  stop('--gwas_group must be specified.\n')
}

# Identify outdir from config file
outdir <- read_param(config = opt$config, 'outdir', return_obj = F)
system(paste0('mkdir -p ', outdir, '/reference/pgs_score_files/', opt$method,'_multi/', opt$gwas_group))

# Create temp directory
tmp_dir<-tempdir()

# Set output prefix
opt$output<-paste0(outdir, '/reference/pgs_score_files/', opt$method,'_multi/', opt$gwas_group,'/ref-', opt$gwas_group)

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'apply_leopard_weights.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

# Read in gwas_groups
gwas_groups <- read_param(config = opt$config, param = 'gwas_groups')
gwas_groups <- gwas_groups[gwas_groups$name == opt$gwas_group,]

# Read in gwas_list
gwas_list <- read_param(config = opt$config, param = 'gwas_list')
gwas_list <- gwas_list[gwas_list$name %in% unlist(strsplit(gwas_groups$gwas, ',')),]

# Split opt$scores
score_files<-paste0(outdir, '/reference/pgs_score_files/', opt$method, '/', gwas_list$name, '/ref-', gwas_list$name, '.score.gz')

#####
# Read in score files and subset pseudo score
#####

log_add(log_file = log_file, message = 'Reading in QuickPRS scores when using full sumstats.')

score_full <- list()
for(i in 1:nrow(gwas_list)){
  param <- find_pseudo(
    config = opt$config,
    gwas = gwas_list$name[i],
    pgs_method = opt$method,
    target_pop = gwas_list$population[i]
  )
  
  score_header <-
    fread(score_files[i], nrows = 1)
  score_cols <-
    which(names(score_header) %in% c('SNP','A1','A2', paste0('SCORE_', param)))
  
  score_full[[gwas_list$population[i]]] <- fread(cmd = 
    paste0(
      "zcat ", score_files[i], " | cut -d' ' -f ", 
      paste0(score_cols, collapse=','))
    )
  
  names(score_full[[gwas_list$population[i]]])[4] <- paste0('SCORE_targ_', gwas_list$population[i])
}

####
# Read in the mixing weights
####

mix_weights <- readRDS(paste0(outdir, '/reference/pgs_score_files/leopard/', opt$gwas_group, '/ref-', opt$gwas_group, '.weights.rds'))

####
# Combine score files
####

log_add(log_file = log_file, message = 'Creating score file.')

# Combine the scores from each population
score_all <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c('SNP','A1','A2'), all = TRUE), score_full)
score_all[is.na(score_all)]<-0

# Read in reference SNP and population data
refdir <- read_param(config = opt$config, param = 'refdir', return_obj = F)
opt$ref_plink_chr <- paste0(refdir, '/ref.chr')
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)
opt$pop_data <- paste0(refdir, '/ref.pop.txt')
pop_data <- read_pop_data(opt$pop_data)
opt$ref_freq_chr <- paste0(refdir, '/freq_files')

# Subset PGS to SNPs in ref (useful when testing)
score_all <- score_all[score_all$SNP %in% ref$SNP,]

# Calculate linear combination of scores using mixing weights for each target population
score_weighted <- score_all
for(targ_pop in gwas_list$population){
  # Read in the .freq file for target population
  freq_data <- read_frq(freq_dir = opt$ref_freq_chr, population = targ_pop, chr = CHROMS)
  
  # Centre SNP-weights for target population
  score_i <- centre_weights(score = score_all, freq = freq_data, ref = ref)
  
  ###
  # Scale weights to give PGS SD of 1 in target population
  ###
  
  # Calculate scores in reference, and scale weights accordingly
  fwrite(score_i, paste0(tmp_dir,'/tmp.',targ_pop,'.score'), col.names=T, sep=' ', quote=F)
  
  # Calc score in target sample
  ref_pgs <- plink_score(pfile = opt$ref_plink_chr, plink2 = 'plink2', keep = pop_data[pop_data$POP == targ_pop, c('FID'), with=F], chr = CHROMS, score = paste0(tmp_dir,'/tmp.',targ_pop,'.score'))
  ref_pgs_scale_i <- score_mean_sd(scores = ref_pgs)
  
  # Rescale SNP-weights according to PGS SD in target
  for(i in gwas_list$population){
    scaling_factor <- 1 / ref_pgs_scale_i$SD[ref_pgs_scale_i$Param == paste0('SCORE_targ_', i)]
    score_i[[paste0('SCORE_targ_', i)]] <- score_i[[paste0('SCORE_targ_', i)]] * scaling_factor
  }
  
  # Linearly combine scores using mixing weights for target population
  score_weighted[[paste0('SCORE_targ_', targ_pop, '_weighted')]] <-
    calculate_weighted_scores(score = score_i,
                              mix_weights = mix_weights,
                              targ_pop = targ_pop)
}

# Only retain weighted score columns
score_weighted <- score_weighted[, c('SNP','A1','A2', names(score_weighted)[grepl('_weighted$', names(score_weighted))]), with=F]

# Reduce number of significant figures to save space
score_weighted[, (4:ncol(score_all)) := lapply(.SD, signif, digits = 7), .SDcols = 4:ncol(score_all)]

fwrite(score_weighted, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

if(file.exists(paste0(opt$output,'.score.gz'))){
  system(paste0('rm ',opt$output,'.score.gz'))
}

system(paste0('gzip ',opt$output,'.score'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()
