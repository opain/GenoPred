#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--ref_plink_chr", action="store", default=NULL, type='character',
		help="Path to per chromosome reference PLINK files [required]"),
make_option("--ref_keep", action="store", default=NULL, type='character',
		help="Keep file to subset individuals in reference for clumping [optional]"),
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
    help="GWAS summary statistics [optional]"),
make_option("--ld_blocks", action="store", default=NULL, type='character',
    help="Path to folder containing LD block information [required]"),
make_option("--dbslmm", action="store", default=NULL, type='character',
    help="Path to DBSLMM directory [required]"),
make_option("--ldsc", action="store", default=NULL, type='character',
    help="Path to LD-score regression binary [required]"),
make_option("--munge_sumstats", action="store", default=NULL, type='character',
    help="Path to munge_sumstats.py script [required]"),
make_option("--ldsc_ref", action="store", default=NULL, type='character',
    help="Path to LD-score regression reference data 'eur_w_ld_chr' [required]"),
make_option("--hm3_snplist", action="store", default=NULL, type='character',
    help="Path to LDSC HapMap3 snplist [required]"),
make_option("--pop_prev", action="store", default=NULL, type='numeric',
    help="Population prevelance (if binary) [optional]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--sample_prev", action="store", default=NULL, type='numeric', 
    help="Sampling ratio in GWAS [optional]")
)

opt = parse_args(OptionParser(option_list = option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../Scripts/functions/misc.R')

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
if(is.null(opt$ld_blocks)){
  stop('--ld_blocks must be specified.\n')
}
if(is.null(opt$dbslmm)){
  stop('--dbslmm must be specified.\n')
}
if(is.null(opt$ldsc)){
  stop('--ldsc must be specified.\n')
}
if(is.null(opt$munge_sumstats)){
  stop('--munge_sumstats must be specified.\n')
}
if(is.null(opt$ldsc_ref)){
  stop('--ldsc_ref must be specified.\n')
}
if(is.null(opt$hm3_snplist)){
  stop('--hm3_snplist must be specified.\n')
}
if(any(!is.null(c(opt$sample_prev, opt$pop_prev))) & any(is.null(c(opt$sample_prev, opt$pop_prev)))){
  stop('If either sample_prev or pop_prev are specified, both must be specified.')
}

# Create output directory
opt$output_dir <- paste0(dirname(opt$output), '/')
system(paste0('mkdir -p ', opt$output_dir))

# Create temp directory
tmp_dir <- tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'dbslmm.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

#####
# Estimate the SNP-heritability using LD-Score Regression
#####

ldsc_h2 <- ldsc(sumstats = opt$sumstats, ldsc = opt$ldsc, munge_sumstats = opt$munge_sumstats, ldsc_ref = opt$ldsc_ref, pop_prev = opt$pop_prev, sample_prev = opt$sample_prev, log_file = log_file)

#####
# Create subset of ref files
#####

if(!is.null(opt$ref_keep)){
  log_add(log_file = log_file, message = 'ref_keep used to subset reference genotype data.')
  plink_subset(bfile = opt$ref_plink_chr, out = paste0(tmp_dir,'/ref_subset_chr'), plink = opt$plink, chr = CHROMS, keep = opt$ref_keep, memory = opt$memory)
  opt$ref_plink_chr_subset<-paste0(tmp_dir,'/ref_subset_chr')
} else {
  opt$ref_plink_chr_subset < - opt$ref_plink_chr
}

#####
# Read in sumstats and insert p-values
#####

log_add(log_file = log_file, message = 'Reading in GWAS and harmonising with reference.')

# Read in, check and format GWAS summary statistics
gwas <- read_sumstats(sumstats = opt$sumstats, chr = CHROMS, log_file = log_file, req_cols = c('CHR','SNP','BP','N','A1','A2','FREQ','BETA','SE','P'))

# Store number of snps and average sample size
nsnp<-nrow(gwas)
gwas_N<-mean(gwas$N)

# Match A1 and A2 match a reference (DBSLMM calls this allele discrepancy)
ref_bim <- read_bim(opt$ref_plink_chr_subset, chr = CHROMS)
gwas <- allele_match(sumstats = gwas, ref_bim = ref_bim, chr = CHROMS)

# Convert to GEMMA format
gwas <- gwas[order(gwas$CHR, gwas$BP),]
gwas$N_MISS <- max(gwas$N) - gwas$N
gwas <- gwas[, c('CHR', 'SNP', 'BP', 'N_MISS', 'N', 'A1', 'A2', 'FREQ', 'BETA', 'SE', 'P'), with=F]
names(gwas)<-c('chr', 'rs', 'ps', 'n_mis', 'n_obs', 'allele1', 'allele0', 'af', 'beta', 'se', 'p_wald')

# Write out formatted sumstats for each chromosome
for(i in CHROMS){
  fwrite(gwas[gwas$chr == i,], paste0(tmp_dir, '/summary_gemma_chr', i, '.assoc.txt'), sep='\t', col.names=F)
}

# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

#####
# Process sumstats using DBSLMM
#####

# This uses the default version of DBSLMM, It appears they now recommend tuning the h2 estimate.
score <- dbslmm(dbslmm = opt$dbslmm, plink = opt$plink, ld_blocks = opt$ld_blocks, chr = CHROMS, bfile = opt$ref_plink_chr_subset, h2 = ldsc_h2, nsnp = nsnp, nindiv = round(gwas_N,0), sumstats = paste0(tmp_dir,'/summary_gemma_chr'), log_file = log_file)

fwrite(score, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

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
ref_pgs <- calc_score(bfile = opt$ref_plink_chr, chr = CHROMS, plink2 = opt$plink2, score = paste0(opt$output,'.score.gz'))

# Calculate scale within each reference population
pop_data <- fread(opt$pop_data)

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

