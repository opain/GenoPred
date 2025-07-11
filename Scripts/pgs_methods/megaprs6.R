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
              help="Path to ldak executable [required]"),
  make_option("--ldak_map", action="store", default=NULL, type='character',
              help="Path to ldak map [required]"),
  make_option("--ldak_tag", action="store", default=NULL, type='character',
              help="Path to ldak tagging data [required]"),
  make_option("--ldak_highld", action="store", default=NULL, type='character',
              help="Path to ldak highld data [required]"),
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
if(is.null(opt$ldak_map)){
  stop('--ldak_map must be specified.\n')
}
if(is.null(opt$ldak_tag)){
  stop('--ldak_tag must be specified.\n')
}
if(is.null(opt$ldak_highld)){
  stop('--ldak_highld must be specified.\n')
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
gwas <- read_sumstats(sumstats = opt$sumstats, chr = CHROMS, log_file = log_file, req_cols = c('CHR','BP','SNP','A1','A2','BETA','SE','N','FREQ'))

# Format for LDAK
snplist <- gwas$SNP
gwas$Z <- gwas$BETA / gwas$SE
gwas$Predictor<-paste0(gwas$CHR, ':', gwas$BP)
gwas<-gwas[,c('Predictor','A1','A2','N','Z','FREQ')]
names(gwas)<-c('Predictor','A1','A2','n','Z','A1Freq')

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
# Format reference for LDAK
############

log_add(log_file = log_file, message = 'Formatting reference for LDAK.')

# Insert CHR:BP IDs
system(paste0("awk < ", tmp_dir, "/ref_merge.bim '{$2=$1\":\"$4;print $0}' > ", tmp_dir, '/tmp.bim; mv ', tmp_dir, '/tmp.bim ', tmp_dir, '/ref_merge.bim'))

# Insert genetic distances
system(paste0(opt$plink, ' --bfile ', tmp_dir, '/ref_merge --cm-map ', opt$ldak_map,'/genetic_map_chr@_combined_b37.txt --make-bed --out ', tmp_dir, '/map'))
system(paste0("cat ", tmp_dir, "/map.bim | awk '{print $2, $3}' > ", tmp_dir, '/map.all'))
system(paste0("awk '(NR==FNR){arr[$1]=$2;next}{print $1, $2, arr[$2], $4, $5, $6}' ", tmp_dir, '/map.all ', tmp_dir, '/ref_merge.bim > ', tmp_dir, '/tmp.bim; mv ', tmp_dir, '/tmp.bim ', tmp_dir, '/ref_merge.bim'))
system(paste0('rm ', tmp_dir, '/map*'))

############
# Estimate Per-Predictor Heritabilities
############
# We will use the BLD-LDAK Model, as recommended for human SNP data

log_add(log_file = log_file, message = 'Estimating per-predictor heritabilities.')

# Calculate LDAK weights
system(paste0(opt$ldak, ' --cut-weights ', tmp_dir,'/sections --bfile ', tmp_dir, '/ref_merge --max-threads ', opt$n_cores))
system(paste0(opt$ldak, ' --calc-weights-all ', tmp_dir,'/sections --bfile ', tmp_dir, '/ref_merge --max-threads ', opt$n_cores))
system(paste0('mkdir ', tmp_dir, '/bld'))
system(paste0('cp ', opt$ldak_tag, '/* ', tmp_dir, '/bld/'))
system(paste0('mv ', tmp_dir, '/sections/weights.short ', tmp_dir,'/bld/bld65'))

# Calculate taggings
if(length(CHROMS) != 1){
  system(paste0(opt$ldak, ' --calc-tagging ', tmp_dir, '/bld.ldak --bfile ', tmp_dir, '/ref_merge --ignore-weights YES --power -.25 --annotation-number 65 --annotation-prefix ', tmp_dir, '/bld/bld --window-cm 1 --save-matrix YES --max-threads ', opt$n_cores))
} else {
  system(paste0(opt$ldak, ' --calc-tagging ', tmp_dir, '/bld.ldak --bfile ', tmp_dir, '/ref_merge --ignore-weights YES --power -.25 --annotation-number 65 --annotation-prefix ', tmp_dir, '/bld/bld --window-cm 1 --chr ', CHROMS, ' --save-matrix YES --max-threads ', opt$n_cores))
}

# Calculate Per-Predictor Heritabilities.
system(paste0(opt$ldak, ' --sum-hers ', tmp_dir, '/bld.ldak --tagfile ', tmp_dir, '/bld.ldak.tagging --summary ', tmp_dir, '/GWAS_sumstats_temp.txt --matrix ', tmp_dir, '/bld.ldak.matrix --max-threads ', opt$n_cores))

ldak_res_her<-fread(paste0(tmp_dir,'/bld.ldak.hers'))

log_add(log_file = log_file, message = paste0('SNP-based heritability estimated to be ',ldak_res_her$Heritability[nrow(ldak_res_her)]," (SD=", ldak_res_her$SD[nrow(ldak_res_her)],")."))

# Identify SNPs in high LD regions
system(paste0(opt$ldak, ' --cut-genes ', tmp_dir, '/highld --bfile ', tmp_dir, '/ref_merge --genefile ', opt$ldak_highld, ' --max-threads ', opt$n_cores))

###################
# Run using full reference.
###################

log_add(log_file = log_file, message = 'Running using full reference.')

# Calculate predictor-predictor correlations
log_add(log_file = log_file, message = 'Calculating predictor-predictor correlations.')
full_cors <- ldak_pred_cor(bfile = paste0(tmp_dir, '/ref_merge'), ldak = opt$ldak, n_cores = opt$n_cores, chr = CHROMS)

# Run MegaPRS
log_add(log_file = log_file, message = paste0('Running MegaPRS: ',opt$prs_model,' model.'))
system(paste0(opt$ldak, ' --mega-prs ', tmp_dir, '/mega_full --model ', opt$prs_model, ' --bfile ', tmp_dir, '/ref_merge --cors ', full_cors, ' --summary ', tmp_dir, '/GWAS_sumstats_temp.txt --allow-ambiguous YES --power -0.25 --max-threads ', opt$n_cores))

# Save the parameters file
system(paste0('cp ', tmp_dir, '/mega_full.parameters ', opt$output, '.model_param.txt'))

# Identify the best fitting model
ldak_res_cors <- fread(paste0(tmp_dir, '/mega_full.cors'), nThread = opt$n_cores)
best_score <- ldak_res_cors[ldak_res_cors$Correlation == max(ldak_res_cors$Correlation),]

# Save the pseudovalidation results
system(paste0('cp ', tmp_dir, '/mega_full.cors ', opt$output, '.pseudoval.txt'))
log_add(log_file = log_file, message = paste0('Model ', gsub('Score_','',best_score$V1[1]),' is identified as the best with correlation of ', best_score$V2))

######
# Format final score file
######

# Read in the scores
score <- fread(paste0(tmp_dir,'/mega_full.effects'), nThread = opt$n_cores)

# Change IDs to RSIDs
ref_pvar <- read_pvar(dat = opt$ref_plink_chr, chr = CHROMS)
ref_pvar$Predictor<-paste0(ref_pvar$CHR,':',ref_pvar$BP)
score<-merge(score, ref_pvar[,c('Predictor','SNP'), with=F], by='Predictor')
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

