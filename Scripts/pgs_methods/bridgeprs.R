#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NULL, type='character',
      help="Path to per chromosome reference PLINK files [required]"),
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

# Split opt$populations
populations<-unlist(strsplit(opt$populations, ','))

#####
# Create LD reference folder
#####
# Bridge requires the LD reference data to be stored in a folder
# Data should be in PLINK1 format, split by chromosome with format "chr<1-22>.ext
# The folder should also include a files called <POP>_ids.txt, which are keep files for each reference population

dir.create(paste0(tmp_dir, '/ref_ld'))
for(i in CHROMS){
  for(j in c('bed','bim','fam')){
    system(paste0(opt$plink2, ' --pfile ', opt$ref_plink_chr, i, ' --make-bed --out ', tmp_dir, '/ref_ld/chr', i))
  }
}

pop_data <- read_pop_data(opt$pop_data)
for(i in unique(pop_data$POP)){
  fwrite(
    pop_data[pop_data$POP == i, c('FID', 'IID'), with = F],
    paste0(tmp_dir, '/ref_ld/', i, '_ids.txt'),
    col.names = F,
    row.names = F,
    quote = F,
    sep = ' '
  )
}

#####
# Prepare sumstats
#####

gwas_N<-NULL
for(i in 1:length(sumstats)){

  log_add(log_file = log_file, message = 'Reading in GWAS.')

  # Read in, check and format GWAS summary statistics
  gwas <- read_sumstats(sumstats = sumstats[i], chr = CHROMS, log_file = log_file, req_cols = c('CHR','SNP','A1','A2','BETA','P','N'))

  # Store average sample size
  gwas_N <- c(gwas_N, round(mean(gwas$N), 0))
  gwas$N<NULL

  # Write sumstats split by chromosome
  for(j in CHROMS){
    tmp<-gwas[gwas$CHR == j,]
    tmp$CHR<-NULL
    fwrite(tmp, paste0(tmp_dir, '/GWAS_sumstats_',i,'_temp_chr', j, '.txt'), sep=' ')
  }

  rm(gwas)
  gc()

}

#####
# Create fake phenotype data
#####

fam<-fread(paste0(tmp_dir, '/ref_ld/chr', CHROMS[1], '.fam'))
names(fam)[1:2]<-c('FID','IID')
pheno<-fam[, 1:2]
pheno$y<-rnorm(nrow(pheno))
write.table(pheno, paste0(tmp_dir, '/fake_pheno.txt'), row.names = F, quote = F)

#####
# Create config files
#####

for(i in 1:length(sumstats)){
  config_i<-c(
    paste0("POP=", populations[i]),
    paste0("LDPOP=", populations[i]),
    paste0("LD_PATH=", tmp_dir, '/ref_ld'),
    paste0("SUMSTATS_PREFIX=", tmp_dir, '/GWAS_sumstats_', i, '_temp_chr'),
    "SUMSTATS_SUFFIX=.txt",
    "SUMSTATS_FIELDS=SNP,A2,A1,P,BETA",
    paste0("SUMSTATS_SIZE=", gwas_N[i]),
    paste0("GENOTYPE_PREFIX=",tmp_dir, '/ref_ld/chr'),
    paste0("PHENOTYPE_FILE=", tmp_dir, '/fake_pheno.txt')
  )

  writeLines(config_i, paste0(tmp_dir, '/', i, '.config'))
}

# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

#####
# Run BridgePRS
#####

system(paste0(opt$bridgeprs_repo, '/bridgePRS pipeline go -o ', tmp_dir,'/out --config_files ', tmp_dir, '/', 1, '.config ', tmp_dir, '/', 2, '.config --phenotype y'))

#####
# Process sumstats using PRS-CSx
#####

# Create a temporary reference bim files for PRS-CSx to match to
pvar <- read_pvar(opt$ref_plink_chr, chr = CHROMS)
pvar$POS<-0
for(i in CHROMS){
  write.table(pvar[pvar$CHR == i, c('CHR','SNP','POS','BP','A1','A2'), with=F], paste0(tmp_dir,'/ref.chr',i,'.bim'), col.names=F, row.names=F, quote=F)
}

rm(pvar)
gc()

# Make a data.frame listing chromosome and phi combinations
jobs<-NULL
for(i in rev(CHROMS)){
  jobs<-rbind(jobs, data.frame(CHR=i, phi=phi_param))
}

# Run using PRScs auto, and specifying a range of global shrinkage parameters
log <- foreach(i = 1:nrow(jobs), .combine = c, .options.multicore = list(preschedule = FALSE)) %dopar% {
  if(jobs$phi[i] == 'auto'){
    system(paste0(opt$prscsx_path, ' --ref_dir=', opt$prscsx_ref_path, '/ --bim_prefix=', tmp_dir,'/ref.chr', jobs$CHR[i], ' --pop=', opt$populations, ' --sst_file=', paste0(paste0(tmp_dir, '/GWAS_sumstats_', 1:length(sumstats),'_temp.txt'), collapse=','),' --n_gwas=', paste(gwas_N, collapse=','), ' --out_dir=', tmp_dir, '/ --out_name=output --chrom=', jobs$CHR[i], ' --meta=True --seed=', opt$seed))
  } else {
    system(paste0(opt$prscsx_path, ' --ref_dir=', opt$prscsx_ref_path, '/ --bim_prefix=', tmp_dir,'/ref.chr', jobs$CHR[i], ' --pop=', opt$populations, ' --phi=', jobs$phi[i], ' --sst_file=', paste0(paste0(tmp_dir, '/GWAS_sumstats_', 1:length(sumstats),'_temp.txt'), collapse=','),' --n_gwas=', paste(gwas_N, collapse=','), ' --out_dir=', tmp_dir, '/ --out_name=output --chrom=', jobs$CHR[i], ' --meta=True --seed=', opt$seed))
  }
}

####
# Combine score files
####

score_all<-NULL
for(pop_i in c(unlist(strsplit(opt$populations, ',')), 'META')){
  score_pop<-NULL
  for(phi_i in phi_param){
    score_phi<-NULL
    for(i in CHROMS){
      score_phi_i<-fread(paste0(tmp_dir,'/output_',pop_i,'_pst_eff_a1_b0.5_phi',phi_i,'_chr',i,'.txt'))
      score_phi<-rbind(score_phi, score_phi_i)
    }
    if(phi_i == phi_param[1]){
      score_phi<-score_phi[,c('V2', 'V4','V5', 'V6'), with=F]
      names(score_phi)<-c('SNP', 'A1', 'A2', paste0('SCORE_',pop_i,'_phi_',phi_i))
    } else {
      score_phi<-score_phi[,'V6', with=F]
      names(score_phi)<-paste0('SCORE_',pop_i,'_phi_',phi_i)
    }
    score_pop<-cbind(score_pop, score_phi)
  }
  if(pop_i == c(unlist(strsplit(opt$populations, ',')), 'META')[1]){
    score_all<-score_pop
  } else {
    score_all<-merge(score_all, score_pop[, !(names(score_pop) %in% c('A1','A2')), with=F], by='SNP', all=T)
  }
}

# Replace NA values with 0
score_all[is.na(score_all)] <- 0

# Flip effects to match reference alleles
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]
score_new <- map_score(ref = ref, score = score_all)

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
