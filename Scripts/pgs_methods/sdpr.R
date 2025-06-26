#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
	make_option("--ref_plink_chr", action="store", default=NULL, type='character',
			help="Path to genome-wide reference PLINK files [required]"),
	make_option("--ref_keep", action="store", default=NULL, type='character',
			help="Keep file to subset individuals in reference for clumping [optional]"),
	make_option("--ref_pcs", action="store", default=NULL, type='character',
	    help="Reference PCs for continuous ancestry correction [optional]"),
	make_option("--gwas_pop", action="store", default=NULL, type='character',
			help="Population of GWAS sample [required]"),
  make_option("--pop_data", action="store", default=NULL, type='character',
      help="File containing the population code and location of the keep file [required]"),
	make_option("--plink2", action="store", default='plink2', type='character',
	    help="Path PLINK v2 software binary [optional]"),
	make_option("--output", action="store", default=NULL, type='character',
			help="Path for output files [required]"),
	make_option("--n_cores", action="store", default=1, type='numeric',
      help="Number of cores to use [optional]"),
	make_option("--pseudo_only", action="store", default=F, type='logical',
      help="Logical indicating whether only pseudovalidated model should be output [optional]"),
	make_option("--sdpr", action="store", default=F, type='character',
      help="Path to SDPR binary [required]"),
	make_option("--test", action="store", default=NA, type='character',
	    help="Specify number of SNPs to include [optional]"),
	make_option("--sumstats", action="store", default=NULL, type='character',
			help="GWAS summary statistics in LDSC format [optional]")
)

opt = parse_args(OptionParser(option_list = option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')
library(lassosum)
library(parallel)
cl <- makeCluster(opt$n_cores)

# Store original working directory
orig_wd<-getwd()

# Check required inputs
if(is.null(opt$ref_plink_chr)){
  stop('--ref_plink_chr must be specified.\n')
}
if(is.null(opt$gwas_pop)){
  stop('--gwas_pop must be specified.\n')
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

# Create output directory
opt$output_dir <- paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'sdpr.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

#####
# Read in sumstats
#####

log_add(log_file = log_file, message = 'Reading in GWAS.')

# Read in, check and format GWAS summary statistics
gwas <- read_sumstats(sumstats = opt$sumstats, chr = CHROMS, log_file = log_file, req_cols = c('CHR','SNP','A1','A2','BETA','SE','N'))

GWAS_CHROMS <- unique(gwas$CHR)
gwas_N <- mean(gwas$N)

# Format for SDPR
gwas$Z <- gwas$BETA/gwas$SE
gwas<-gwas[, c('SNP','A1','A2','Z'), with = F]
fwrite(gwas, paste0(tmp_dir,'/sumstats.txt'), row.names=F, quote=F, sep=' ', na='NA')
write.table(gwas$SNP, paste0(tmp_dir,'/sumstats.extract'), col.names=F, row.names=F, quote=F)

###
# Run SDPR
###

dir.create(paste0(tmp_dir, '/ref_ld'))
dir.create(paste0(tmp_dir, '/ref_mat'))
dir.create(paste0(tmp_dir, '/result'))

log_add(log_file = log_file, message = 'Running SDPR.')

score <- NULL
h2_total <- 0

run_sdpr <- function(chr, r2 = 0.1) {
  # Build LD reference
  system(paste0(
    opt$sdpr, ' -make_ref -ref_prefix ', tmp_dir, '/ref_ld/chr', chr,
    ' -r2 ', r2,
    ' -chr ', chr,
    ' -ref_dir ', tmp_dir, '/ref_mat/'
  ))
  
  # Run SDPR and return output
  system(paste0(
    opt$sdpr,
    ' -mcmc -ref_dir ', tmp_dir, '/ref_mat/',
    ' -ss ', tmp_dir, '/sumstats.txt',
    ' -N ', round(gwas_N, 0),
    ' -chr ', chr,
    ' -out ', tmp_dir, '/result/SDPR_chr', chr, '.txt'
  ), intern = TRUE)
}

for (i in GWAS_CHROMS) {
  message("Running SDPR on chr", i)
  
  # Build LD reference
  system(paste0(
    opt$plink2, ' --pfile ', opt$ref_plink_chr, i,
    ' --keep ', opt$ref_keep,
    ' --extract ', tmp_dir, '/sumstats.extract',
    ' --make-bed --out ', tmp_dir, '/ref_ld/chr', i
  ))

  r2<-0.1
  while(i){
    sdpr_output <- run_sdpr(i, r2 = r2)
    h2_line <- sdpr_output[grepl("^h2:", sdpr_output)]
    if (!grepl('h2: -nan', h2_line)) {
      break
    }
    
    log_add(log_file = log_file, message = paste0("SDPR failed on chr ", i, " even with r2=", r2, "."))
    r2 <- r2 + 0.1
  }
  
  # Parse and sum heritability
  h2_val <- as.numeric(sub("h2: ([0-9.]+).*", "\\1", h2_line))
  h2_total <- h2_total + h2_val
  
  # Load the score output
  score_file <- paste0(tmp_dir, '/result/SDPR_chr', i, '.txt')
  score <- rbind(score, fread(score_file))

}

log_add(log_file = log_file, message = paste0('SNP-based heritability = ', h2_total, '.'))

# Insert A2 info
ref <- read_pvar(opt$ref_plink_chr, chr = CHROMS)[, c('SNP','A1','A2'), with=F]
score <- merge(score, ref, by = c('SNP','A1'))
names(score)[names(score) == 'beta'] <- 'SCORE_SDPR'
score <- score[, c('SNP','A1','A2','SCORE_SDPR'), with=F]

# Flip effects to match reference alleles
score <- map_score(ref = ref, score = score)

# Reduce number of significant figures to save space
score[, (4:ncol(score)) := lapply(.SD, signif, digits = 7), .SDcols = 4:ncol(score)]

fwrite(score, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

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
