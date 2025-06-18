#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NULL, type='character',
              help="Path to per chromosome reference PLINK files [required]"),
  make_option("--ref_keep", action="store", default=NULL, type='character',
              help="Path to keep file for reference [optional]"),
  make_option("--ref_pcs", action="store", default=NULL, type='character',
              help="Reference PCs for continuous ancestry correction [optional]"),
  make_option("--pop_data", action="store", default=NULL, type='character',
              help="File containing the population code and location of the keep file [required]"),
  make_option("--plink", action="store", default='plink', type='character',
              help="Path PLINK v1.9 software binary [optional]"),
  make_option("--plink2", action="store", default='plink2', type='character',
              help="Path PLINK v2 software binary [optional]"),
  make_option("--output", action="store", default=NULL, type='character',
              help="Path for output files [required]"),
  make_option("--memory", action="store", default=5000, type='numeric',
              help="Memory limit [optional]"),
  make_option("--sumstats", action="store", default=NULL, type='character',
              help="GWAS summary statistics [required]"),
  make_option("--ldak", action="store", default=NULL, type='character',
              help="Path to ldak v6.1 executable [required]"),
  make_option("--breakpoints", action="store", default=NULL, type='character',
              help="Path to breakpoints file [required]"),
  make_option("--quickprs_ldref", action="store", default=NA, type='character',
              help="Path to folder containing ldak quickprs reference [required]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
              help="Number of cores for parallel computing [optional]"),
  make_option("--prs_model", action="store", default='mega', type='character',
              help="Model used for deriving SNP-weights [optional]"),
  make_option("--test", action="store", default=NA, type='character',
              help="Specify number of SNPs to include [optional]")
)

opt = parse_args(OptionParser(option_list = option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')

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
if(is.null(opt$ldak)){
  stop('--ldak must be specified.\n')
}

if(!is.na(opt$quickprs_ldref) && opt$quickprs_ldref == 'NA'){
  opt$quickprs_ldref<-NA
}

if(is.na(opt$quickprs_ldref)){
  if(is.null(opt$breakpoints)){
    stop('--breakpoints must be specified when --quickprs_ldref is not specified.\n')
  }
}

# Create output directory
opt$output_dir <- paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'megaprs.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

#####
# Format the sumstats
#####

log_add(log_file = log_file, message = 'Reading in GWAS.')

# Read in, check and format GWAS summary statistics
gwas <- read_sumstats(sumstats = opt$sumstats, chr = CHROMS, log_file = log_file, req_cols = c('CHR','BP','SNP','A1','A2','BETA','SE','N','FREQ','REF.FREQ'))

# Format for LDAK
snplist <- gwas$SNP
gwas$Z <- gwas$BETA / gwas$SE
gwas<-gwas[,c('SNP','A1','A2','N','Z')]
names(gwas)<-c('Predictor','A1','A2','n','Z')

if(!is.na(opt$quickprs_ldref)){
  
  log_add(log_file = log_file, message = 'Using precomputed predictor-predictor correlations.')
  
  # Check overlap between GWAS and LDAK reference
  ldak_hm3_file <- list.files(opt$quickprs_ldref)
  ldak_hm3_file <- ldak_hm3_file[grepl('.cors.bim', ldak_hm3_file)][1]
  ldak_hm3 <- fread(paste0(opt$quickprs_ldref, '/', ldak_hm3_file))
  ldak_hm3 <- ldak_hm3[ldak_hm3$V1 %in% CHROMS,]
  ref_overlap <- sum(gwas$Predictor %in% ldak_hm3$V2) / nrow(ldak_hm3)
  
  log_add(log_file = log_file, message = paste0('GWAS-reference overlap is ', round(ref_overlap * 100, 2), '%.'))

  # Subset GWAS to LDAK reference data
  gwas <- gwas[gwas$Predictor %in% ldak_hm3$V2, ]

}

fwrite(gwas, paste0(tmp_dir,'/GWAS_sumstats_temp.txt'), sep=' ')

###
# Merge the per chromosome reference genetic data and subset opt$ref_keep
###

log_add(log_file = log_file, message = 'Merging per chromosome reference data.')

# Save in plink1 format for MegaPRS
plink_merge(pfile = opt$ref_plink_chr, chr = CHROMS, plink2 = opt$plink2, keep = opt$ref_keep, extract = snplist, make_bed =T, out = paste0(tmp_dir, '/ref_merge'))

# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

############
# Create LD matrices
############

if(is.na(opt$quickprs_ldref)){
  
  log_add(log_file = log_file, message = 'Calculating predictor-predictor correlations.')
  
  system(paste0(
    opt$ldak, ' ',
    '--calc-cors ', tmp_dir, '/full_cors ',
    '--bfile ', tmp_dir, '/ref_merge ', 
    '--break-points ', opt$breakpoints, ' ',
    '--max-threads ', opt$n_cores))
  
  ld_ref_path <- paste0(tmp_dir, '/full_cors')

} else {

  ref_files<-list.files(opt$quickprs_ldref)
  cor_file_prefix<-gsub('.cors.bin','',ref_files[grepl('.cors.bin',ref_files) & !grepl('subset', ref_files)])
  ld_ref_path <- paste0(opt$quickprs_ldref,'/',cor_file_prefix)
  
}

############
# Run MegaPRS
############

# Note. MegaPRS by default runs pseudo CV to identify the best model, and then 
# only produces that model using the full sumstats. I would like to output 
# weights from all models, but also use pseudo CV to identify the best model.
# As a temporary work around I will run MegaPRS twice with and without the CV
# --skip-cv YES

log_add(log_file = log_file, message = paste0('Running MegaPRS: ',opt$prs_model,' model.'))

system(paste0(
  opt$ldak, ' ',
  '--mega-prs ', tmp_dir, '/mega_full ',
  '--model ', opt$prs_model, ' ',
  '--cors ', ld_ref_path, ' ',
  '--summary ', tmp_dir, '/GWAS_sumstats_temp.txt ',
  '--skip-cv NO ',
  '--check-sums NO ',
  '--max-threads ', opt$n_cores))

system(paste0(
  opt$ldak, ' ',
  '--mega-prs ', tmp_dir, '/mega_full ',
  '--model ', opt$prs_model, ' ',
  '--cors ', ld_ref_path, ' ',
  '--summary ', tmp_dir, '/GWAS_sumstats_temp.txt ',
  '--skip-cv YES ',
  '--check-sums NO ',
  '--max-threads ', opt$n_cores))

# Save the parameters file
system(paste0('cp ', tmp_dir, '/mega_full.parameters ', opt$output, '.model_param.txt'))

# Save the pseudosummary results
system(paste0('cp ',tmp_dir,'/mega_full.cors ',opt$output,'.pseudoval.txt'))

# Identify the best fitting model
ldak_res_cors <- fread(paste0(tmp_dir, '/mega_full.cors'), nThread = opt$n_cores)
best_score <- ldak_res_cors[which.max(ldak_res_cors$Correlation),]

log_add(log_file = log_file, message = paste0('Model ', gsub('Score_','',best_score$Model[1]),' is identified as the best with correlation of ', best_score$Correlation[1]))

######
# Format final score file
######

# Read in the scores
score <- fread(paste0(tmp_dir,'/mega_full.effects'), nThread = opt$n_cores)
names(score)[1]<-'SNP'
score<-score[, c('SNP', 'A1', 'A2', names(score)[grepl('Model', names(score))]), with=F]
names(score)[grepl('Model', names(score))]<-paste0('SCORE_ldak_',names(score)[grepl('Model', names(score))])

# Flip effects to match reference alleles
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]
score_new <- map_score(ref = ref, score = score)

fwrite(score_new, paste0(opt$output, '.score'), col.names=T, sep=' ', quote=F)

if(file.exists(paste0(opt$output,'.score.gz'))){
  system(paste0('rm ',opt$output,'.score.gz'))
}

system(paste0('gzip ',opt$output,'.score'))

# Record end time of test
if(!is.na(opt$test)){
  test_finish(log_file = log_file, test_start.time = test_start.time)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()

