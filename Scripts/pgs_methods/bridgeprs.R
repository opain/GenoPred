#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--ref_plink_chr", action="store", default=NA, type='character',
      help="Path to per chromosome reference PLINK files [required]"),
  make_option("--ref_pcs", action="store", default=NULL, type='character',
      help="Reference PCs for continuous ancestry correction [optional]"),
  make_option("--pop_data", action="store", default=NULL, type='character',
      help="File containing the population code and location of the keep file [required]"),
  make_option("--plink2", action="store", default='plink2', type='character',
      help="Path PLINK v2 software binary [optional]"),
  make_option("--bridgeprs_repo", action="store", default=NULL, type='character',
      help="BridgePRS repo path [optional]"),
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
if(is.null(opt$bridgeprs_repo)){
  stop('--bridgeprs_repo must be specified.\n')
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
pop_data$FID <- 0
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
  
  # FOR testing
  gwas <- gwas[sample(1:nrow(gwas), 10000),]

  # Write sumstats split by chromosome
  for(j in CHROMS){
    tmp<-gwas[gwas$CHR == j,]
    tmp$CHR<-NULL
    fwrite(tmp, paste0(tmp_dir, '/GWAS_sumstats_',populations[i],'_temp_chr', j, '.txt'), sep=' ')
  }
  
  # Write QC snplist
  fwrite(gwas[,'SNP'], paste0(tmp_dir, '/GWAS_sumstats_',populations[i],'_temp.qc.snplist'), col.names = F, quote = F)

  rm(gwas)
  gc()

}

#####
# Run BridgePRS
#####

# set fst
if(any(populations == 'EUR') & any(populations == 'AFR')){
  fst <- 0.15
}
if(any(populations == 'EUR') & any(populations == 'EAS')){
  fst <- 0.11
}

system(paste0('rm -r ', tmp_dir,'/bridge_out'))
dir.create(paste0(tmp_dir,'/bridge_out'))

system(paste0(
  opt$bridgeprs_repo, '/src/Bash/BridgePRS_sumstat.sh ', opt$bridgeprs_repo, '/src/Rscripts ',
  '--strand_check 1 ',
  '--outdir ', tmp_dir,'/bridge_out ',
  '--by_chr 1 ',
  '--by_chr_sumstats .txt ',
  '--pop1 ', populations[1], ' ',
  '--pop2 ', populations[2], ' ',
  '--fst ', fst, ' ',
  '--N_pop1 ', gwas_N[1], ' ',
  '--N_pop2 ', gwas_N[2], ' ',
  '--pop1_sumstats ', tmp_dir, '/GWAS_sumstats_', populations[1], '_temp_chr ',
  '--pop1_qc_snplist ', tmp_dir, '/GWAS_sumstats_',populations[1],'_temp.qc.snplist ',
  '--pop2_sumstats ', tmp_dir, '/GWAS_sumstats_', populations[2], '_temp_chr ',
  '--pop2_qc_snplist ', tmp_dir, '/GWAS_sumstats_',populations[2],'_temp.qc.snplist ',
  '--pop1_ld_bfile ', tmp_dir, '/ref_ld/chr ',
  '--pop1_ld_ids ', tmp_dir, '/ref_ld/', populations[1], '_ids.txt ',
  '--pop2_ld_bfile ', tmp_dir, '/ref_ld/chr ',
  '--pop2_ld_ids ', tmp_dir, '/ref_ld/', populations[2], '_ids.txt ',
  '--sumstats_snpID SNP ',
  '--sumstats_p P ',
  '--sumstats_beta BETA ',
  '--sumstats_allele1 A1 ',
  '--sumstats_allele0 A2 ',
  '--pheno_name y ',
  '--n_cores ', opt$n_cores, ' ',
  '--prop_train 0.6 ',
  '--prop_test 0.3 ',
  '--do_block_pop1 1 ',
  '--do_block_pop2 1 ',
  '--do_sumstat_pop1 1 ',
  '--do_sumstat_pop2 1 ',
  '--do_clump_pop1 1 ',
  '--do_est_beta_pop1 1 ',
  '--do_sumstat_ensembl_pop1 1 ',
  '--do_est_beta_pop1_precision 1 ',
  '--do_est_beta_InformPrior 1 ',
  '--do_clump_pop2 1 ',
  '--do_est_beta_pop2 1 ',
  '--do_sumstat_ensembl_pop2 1 ',
  '--n_folds 10 ',
  ' > ',tmp_dir,'/bridge_log.txt 2>&1'
))

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

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()
