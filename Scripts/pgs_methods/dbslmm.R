#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--ref_plink_chr", action="store", default=NULL, type='character',
		help="Path to per chromosome reference PLINK files [required]"),
make_option("--ref_keep", action="store", default=NULL, type='character',
    help="Keep file to subset individuals in reference for clumping [optional]"),
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
    help="GWAS summary statistics [optional]"),
make_option("--ld_blocks", action="store", default=NULL, type='character',
    help="Path to folder containing LD block information [required]"),
make_option("--dbslmm", action="store", default=NULL, type='character',
    help="Path to DBSLMM directory [required]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--n_cores", action="store", default=1, type='numeric',
    help="Number of cores for parallel computing [optional]"),
make_option("--h2f", action="store", default='0.8,1,1.2', type='character',
    help="Folds of SNP-based heritability [optional]"),
make_option("--h2", action="store", default=NULL, type='numeric',
    help="SNP-based heritability [optional]")

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
if(is.null(opt$ld_blocks)){
  stop('--ld_blocks must be specified.\n')
}
if(is.null(opt$dbslmm)){
  stop('--dbslmm must be specified.\n')
}
if(is.null(opt$h2)){
  stop('--h2 must be specified.\n')
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

# Format the h2f parameter
opt$h2f <- as.numeric(unlist(strsplit(opt$h2f, ',')))

if(file.exists(opt$h2)){
  opt$h2<-read.table(opt$h2, header=F)$V1
} else {
  opt$h2<-0.05
}

if(opt$h2 > 1){
  opt$h2 <- 1
  log_add(log_file = log_file, message = 'SNP-h2 was set to 1.')
}

if(opt$h2 < 0.05){
  opt$h2 <- 0.05
  log_add(log_file = log_file, message = 'SNP-h2 was set to 0.05.')
}

if(any(opt$h2*opt$h2f > 1)){
  opt$h2 <- opt$h2*(1/max(opt$h2f*opt$h2))
  log_add(log_file = log_file, message = paste0('SNP-h2 was set to ',opt$h2,' to avoid SNP-h2*h2f > 1.'))
}

# Read in and write out out sumstats to avoid munge error
gwas <- read_sumstats(sumstats = opt$sumstats, chr = CHROMS, log_file = log_file, req_cols = c('CHR','SNP','BP','N','A1','A2','FREQ','BETA','SE','P'))
fwrite(gwas, paste0(tmp_dir,'/sumstats.gz'), row.names=F, quote=F, sep=' ', na='NA')
opt$sumstats<-paste0(tmp_dir,'/sumstats.gz')
GWAS_CHROMS <- unique(gwas$CHR)

#####
# Create subset of ref files
#####
# Save in plink1 format for DBSLMM
if(!is.null(opt$ref_keep)){
  log_add(log_file = log_file, message = 'ref_keep used to subset reference genotype data.')
  plink_subset(pfile = opt$ref_plink_chr, make_bed = T, out = paste0(tmp_dir,'/ref_subset_chr'), extract = gwas$SNP, plink2 = opt$plink2, chr = GWAS_CHROMS, keep = opt$ref_keep, memory = opt$memory)
  opt$ref_plink_chr_subset<-paste0(tmp_dir,'/ref_subset_chr')
} else {
  plink_subset(pfile = opt$ref_plink_chr, make_bed = T, out = paste0(tmp_dir,'/ref_subset_chr'), extract = gwas$SNP, plink2 = opt$plink2, chr = GWAS_CHROMS, memory = opt$memory)
  opt$ref_plink_chr_subset<-paste0(tmp_dir,'/ref_subset_chr')
}

#####
# Read in sumstats and insert p-values
#####

log_add(log_file = log_file, message = 'Reading in GWAS and harmonising with reference.')

# Store number of snps and average sample size
nsnp<-nrow(gwas)
gwas_N<-mean(gwas$N)

# Match A1 and A2 match a reference (DBSLMM calls this allele discrepancy)
ref_bim <- read_bim(opt$ref_plink_chr_subset, chr = GWAS_CHROMS)
gwas <- allele_match(sumstats = gwas, ref_bim = ref_bim, chr = GWAS_CHROMS)

# Convert to GEMMA format
gwas <- gwas[order(gwas$CHR, gwas$BP),]
gwas$N_MISS <- max(gwas$N) - gwas$N
gwas <- gwas[, c('CHR', 'SNP', 'BP', 'N_MISS', 'N', 'A1', 'A2', 'FREQ', 'BETA', 'SE', 'P'), with=F]
names(gwas)<-c('chr', 'rs', 'ps', 'n_mis', 'n_obs', 'allele1', 'allele0', 'af', 'beta', 'se', 'p_wald')

# Write out formatted sumstats for each chromosome
for(i in GWAS_CHROMS){
  fwrite(gwas[gwas$chr == i,], paste0(tmp_dir, '/summary_gemma_chr', i, '.assoc.txt'), sep='\t', col.names=F)
}

# Record start time for test
if(!is.na(opt$test)){
  test_start.time <- test_start(log_file = log_file)
}

#####
# Process sumstats using DBSLMM
#####

# We will use the tuning version of DBSLMM, where h2 is multiplied by a factor of 0.8, 1, and 1.2
score <-
  dbslmm(
    dbslmm = opt$dbslmm,
    plink = opt$plink,
    ld_blocks = opt$ld_blocks,
    chr = GWAS_CHROMS,
    bfile = opt$ref_plink_chr_subset,
    h2 = opt$h2,
    h2f = opt$h2f,
    nsnp = nsnp,
    ncores = opt$n_cores,
    nindiv = round(gwas_N, 0),
    sumstats = paste0(tmp_dir, '/summary_gemma_chr'),
    log_file = log_file
  )

# Flip effects to match reference alleles
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]
score_new <- map_score(ref = ref, score = score)

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

