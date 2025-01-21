#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--config", action="store", default=NULL, type='character',
      help="Pipeline configuration file. Required when pseudo_only is TRUE [optional]"),
  make_option("--ref_plink_chr", action="store", default=NA, type='character',
      help="Path to per chromosome reference PLINK files [required]"),
  make_option("--ref_pcs", action="store", default=NULL, type='character',
      help="Reference PCs for continuous ancestry correction [optional]"),
  make_option("--pop_data", action="store", default=NULL, type='character',
      help="File containing the population code and location of the keep file [required]"),
  make_option("--ref_keep_dir", action="store", default=NULL, type='character',
      help="Directory continaing reference keep files [required]"),
  make_option("--plink1", action="store", default='plink', type='character',
      help="Path PLINK v1.9 software binary [optional]"),
  make_option("--plink2", action="store", default='plink2', type='character',
      help="Path PLINK v2 software binary [optional]"),
  make_option("--output", action="store", default=NULL, type='character',
      help="Path for output files [required]"),
  make_option("--memory", action="store", default=5000, type='numeric',
      help="Memory limit [optional]"),
  make_option("--sumstats", action="store", default=NULL, type='character',
      help="Comma-seperated list of GWAS summary statistics [required]"),
  make_option("--scores", action="store", default=NULL, type='character',
      help="Comma-seperated list of score files [required]"),
  make_option("--populations", action="store", default=NULL, type='character',
      help="Comma-seperated list of population codes matching GWAS [required]"),
  make_option("--retain_nonoverlapping", action="store", default=T, type='character',
      help="Logical indicating whether or not to retain the original BETA if variant is missing in target GWAS [required]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
      help="Number of cores for parallel computing [optional]"),
  make_option("--test", action="store", default=NA, type='character',
      help="Specify number of SNPs to include [optional]"),
  make_option("--pseudo_only", action="store", default=T, type='character',
      help="Apply TLPRS to model selected by pseudovalidation only [optional]"),
  make_option("--seed", action="store", default=1, type='numeric',
      help="Seed number for PRScs [optional]")
)

opt = parse_args(OptionParser(option_list = option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')
library(lassosum)
library(TLPRS)

library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

# Check required inputs
if(is.null(opt$config) && opt$pseudo_only){
  stop('--config must be specified when --pseudo_only is TRUE.\n')
}
if(is.null(opt$ref_plink_chr)){
  stop('--ref_plink_chr must be specified.\n')
}
if(is.null(opt$sumstats)){
  stop('--sumstats must be specified.\n')
}
if(is.null(opt$scores)){
  stop('--sumstats must be specified.\n')
}
if(is.null(opt$pop_data)){
  stop('--pop_data must be specified.\n')
}
if(is.null(opt$ref_keep_dir)){
  stop('--ref_keep_dir must be specified.\n')
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

# Split opt$scores
scores<-unlist(strsplit(opt$scores, ','))
log_add(log_file = log_file, message = paste0(length(scores), ' sets of scores have been provided.'))

# Split opt$populations
populations<-unlist(strsplit(opt$populations, ','))

######
# Merge reference data
######

for(i in 1:length(populations)) {
  plink_merge(
    pfile = opt$ref_plink_chr,
    chr = CHROMS,
    plink2 = opt$plink2,
    keep = paste0(opt$ref_keep_dir, '/', populations[i], '.keep'),
    make_bed = T,
    out = paste0(tmp_dir, '/', populations[i], '_ref_merge')
  )
}

# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

######
# Run TL-PRS
######
# We are going to use code within the TL_PRS function to generate the BETAs across across a range of gradients, avoid the validation step

# Read in reference SNP data for harmonising across GWAS
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]

# Run for each parameter in each score file
tl_betas_list<-list()
for(i in 1:length(populations)){

  # Read in target sumstats
  target_gwas <- read_sumstats(sumstats = sumstats[i], chr = CHROMS, log_file = log_file, req_cols = c('SNP', 'A1', 'BETA', 'P', 'N'))
  names(target_gwas)[names(target_gwas) == 'BETA']<-'beta'
  names(target_gwas)[names(target_gwas) == 'P']<-'p'

  # Identify LD block data to be used
  if(populations[i] %in% c('EUR','AFR')){
    ld_block_dat <- paste0(populations[i],'.hg19')
  }
  if(populations[i] %in% 'EAS'){
    ld_block_dat <- 'ASN.hg19'
  }
  if(populations[i] %in% c('AMR','SAS')){
    ld_block_dat <- 'EUR.hg19'
    log_add(log_file = log_file, message = 'Using LD block data for EUR.')
  }

  # Read in PGS score file (set up for pairs of GWAS/populations only)
  score_file <- fread(scores[-i])

  if(opt$pseudo_only){
    method <- sub('/.*','',gsub('.*pgs_score_files/','', scores[-i]))
    gwas <- sub('/.*','',gsub(paste0('.*/',method,'/'),'', scores[-i]))
    param <- find_pseudo(
      config = opt$config,
      gwas = gwas,
      pgs_method = method,
      target_pop = populations[i]
    )
    score_file <- score_file[, c('SNP', 'A1', 'A2', paste0('SCORE_', param)), with = F]
    log_add(log_file = log_file, message = 'Using pseudovalidated PGS only.')
  }

  names(score_file)<-gsub('^SCORE','Beta', names(score_file))

  target_gwas_j=merge(score_file, target_gwas, by="SNP",sort=F)

  if (sum(target_gwas_j$p <= 1E-320)>0) {
    target_gwas_j$p[target_gwas_j$p <= 1E-320] = 1E-320
  }

  target_gwas_j$cor = lassosum::p2cor(
    p = target_gwas_j$p,
    n = median(target_gwas_j$N, na.rm = T),
    sign = target_gwas_j$beta
  )

  flag = which(target_gwas_j$A1.x != target_gwas_j$A1.y)
  if (length(flag) > 0) {
    target_gwas_j$cor[flag] = -target_gwas_j$cor[flag]
  }
  target_gwas_j=target_gwas_j[,c("SNP","A1.x", names(target_gwas_j)[grepl('^Beta', names(target_gwas_j))], "cor"), with=F]
  colnames(target_gwas_j)[2]="A1"
  gc()

  beta_list = as.data.frame(
    TLPRS:::PRStr_calculation2(
      sum_stats_target = target_gwas_j,
      train_file = paste0(tmp_dir, '/', populations[i], '_ref_merge'),
      sum_stats = score_file,
      LDblocks = ld_block_dat,
      plink=opt$plink1,
      temp.file = paste0(tmp_dir, '/', populations[i], '_step1')
    )
  )

  # Flip effects to correspond to original A1
  beta.info<-beta_list[, 1:9]
  beta_list<-beta_list[, -1:-9]
  flip<-beta.info$V5 != beta.info$A1
  beta_list[flip,]<- -beta_list[flip,]

  for (k in 1:ncol(beta_list)){
    sdtemp=sd(beta_list[,k],na.rm=T)
    if (sdtemp>1){
      beta_list[,k:ncol(beta_list)]=1
    }
  }

  beta_list=beta_list/beta.info$sd

  colnames(beta.info)[1:2]=c("SNP","A1")
  beta.info<-beta.info[, 1:2]

  beta_list<-data.table(cbind(beta.info, beta_list))
  beta_list<-merge(score_file[, c('SNP','A1','A2'), with=F], beta_list, by=c('SNP','A1'))

  names(beta_list)<-gsub('^Beta', paste0('SCORE_targ_', populations[i]), names(beta_list))

  if(opt$retain_nonoverlapping){
    # Insert original BETA if variant not present in target GWAS
    miss_snps <- score_file[!(score_file$SNP %in% beta_list$SNP),]
    beta_list <- merge(beta_list, miss_snps, by = c('SNP','A1','A2'), all=T)
    for(j in gsub('Beta_', '', names(score_file)[-1:-3])){
      beta_list[!is.na(get(paste0('Beta_', j))),
                (which(grepl('SCORE', names(beta_list)) & grepl(j, names(beta_list))))] <-
        beta_list[[paste0('Beta_', j)]][!is.na(beta_list[[paste0('Beta_', j)]])]
    }
    beta_list<-beta_list[, !grepl('Beta_', names(beta_list)), with=F]
  }

  # Flip effects to match reference alleles
  beta_list <- map_score(ref = ref, score = beta_list)

  tl_betas_list[[i]]<-beta_list

}

tl_betas_all<-Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c('SNP','A1','A2'), all = TRUE, sort = F), tl_betas_list)

# Reduce number of significant figures to save space
tl_betas_all[, (4:ncol(tl_betas_all)) := lapply(.SD, signif, digits = 7), .SDcols = 4:ncol(tl_betas_all)]

fwrite(tl_betas_all, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

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
