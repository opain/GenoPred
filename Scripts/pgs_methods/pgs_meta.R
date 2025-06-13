#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--output", action="store", default='NA', type='character',
              help="Path for output files [required]"),
  make_option("--plink2", action="store", default='plink2', type='character',
              help="Path PLINK v2 software binary [optional]"),
  make_option("--ref_plink_chr", action="store", default=NA, type='character',
              help="Path to per chromosome reference PLINK files [required]"),
  make_option("--pop_data", action="store", default=NULL, type='character',
              help="File containing the population code and location of the keep file [required]"),
  make_option("--method", action="store", default=NULL, type='character',
              help="PGS method [required]"),
  make_option("--sumstats", action="store", default=NULL, type='character',
              help="Comma-seperated list of GWAS summary statistics [required]"),
  make_option("--scores", action="store", default=NULL, type='character',
              help="Comma-seperated list of score files [required]"),
  make_option("--populations", action="store", default=NULL, type='character',
              help="Comma-seperated list of population codes matching GWAS [required]"),
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
if(is.null(opt$pop_data)){
  stop('--pop_data must be specified.\n')
}
if(is.null(opt$sumstats)){
  stop('--sumstats must be specified.\n')
}
if(is.null(opt$method)){
  stop('--method must be specified.\n')
}
if(is.null(opt$scores)){
  stop('--scores must be specified.\n')
}
if(is.null(opt$output)){
  stop('--output must be specified.\n')
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
log_header(log_file = log_file, opt = opt, script = 'pgsmeta', start.time = start.time)

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

# Split opt$scores
score_files<-unlist(strsplit(opt$scores, ','))

###
# Read in GWAS sumstats
###

sumstats_list<-list()
for(i in 1:length(sumstats)){
  sumstats_list[[i]] <- read_sumstats(sumstats = sumstats[i], chr = CHROMS, log_file = log_file, req_cols = c('SNP','N'))
}

###
# Read in score files
###

score_list<-list()
for(i in 1:length(score_files)){
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
  
  score_list[[i]] <- fread(score_files[i])
}

###
# Merge sumstats and respective scores
###

both_list<-list()
for(i in 1:length(sumstats)){
  both_list[[i]] <- merge(score_list[[i]], sumstats_list[[i]], by = c('SNP'), sort = F, all.x = T)
}

###
# Meta analyse SNP effects
###

# Stopped here as realised we need to extract pseudoval PGS requiring config
meta_analyse_snp_weights <- function(both_list, score_col = "SCORE_quickprs", weight_col = "N") {
  # Merge all data frames by SNP (full outer join)
  merged <- Reduce(function(x, y) merge(x, y, by = "SNP", all = TRUE), score_list)
  
  k <- length(both_list)
  score_cols <- paste0(score_col, "_", 1:k)
  weight_cols <- paste0(weight_col, "_", 1:k)
  
  # Ensure numeric
  for (i in 1:k) {
    merged[[score_cols[i]]] <- as.numeric(merged[[score_cols[i]]])
    merged[[weight_cols[i]]] <- as.numeric(merged[[weight_cols[i]]])
  }
  
  # Meta-analyse each SNP using available weights
  numerator <- rowSums(sapply(1:k, function(i) {
    score <- merged[[score_cols[i]]]
    weight <- merged[[weight_cols[i]]]
    score * weight * !is.na(score) * !is.na(weight)
  }), na.rm = TRUE)
  
  denominator <- rowSums(sapply(1:k, function(i) {
    weight <- merged[[weight_cols[i]]]
    !is.na(merged[[score_cols[i]]]) * !is.na(weight) * weight
  }), na.rm = TRUE)
  
  # Compute meta score only where denominator > 0
  merged$meta_score <- ifelse(denominator > 0, numerator / denominator, NA)
  
  return(merged[, c("SNP", "meta_score")])
}

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

# Flip effects to match reference alleles
score_all <- map_score(ref = ref, score = score_all)

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
score_weighted[, (4:ncol(score_weighted)) := lapply(.SD, signif, digits = 7), .SDcols = 4:ncol(score_weighted)]

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
