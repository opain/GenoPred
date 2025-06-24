#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--ref_plink_chr", action="store", default=NULL, type='character',
    help="Path to per chromosome reference PLINK files [required]"),
make_option("--output", action="store", default=NULL, type='character',
		help="Path for output files [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]"),
make_option("--sumstats", action="store", default=NULL, type='character',
    help="GWAS summary statistics [optional]"),
make_option("--ldsc", action="store", default=NULL, type='character',
    help="Path to LD-score regression binary [required]"),
make_option("--munge_sumstats", action="store", default=NULL, type='character',
    help="Path to munge_sumstats.py script [required]"),
make_option("--ld_scores", action="store", default=NULL, type='character',
    help="Path to genome-wide ld scores [required]"),
make_option("--hm3_snplist", action="store", default=NULL, type='character',
    help="Path to LDSC HapMap3 snplist [required]"),
make_option("--hm3_no_mhc", action="store", default=F, type='logical',
    help="Logical indicating whether MHC region should be removed for LDSC analysis [required]"),
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
source('../functions/misc.R')
source_all('../functions')

# Check required inputs
if(is.null(opt$ref_plink_chr)){
  stop('--ref_plink_chr must be specified.\n')
}
if(is.null(opt$sumstats)){
  stop('--sumstats must be specified.\n')
}
if(is.null(opt$output)){
  stop('--output must be specified.\n')
}
if(is.null(opt$ldsc)){
  stop('--ldsc must be specified.\n')
}
if(is.null(opt$munge_sumstats)){
  stop('--munge_sumstats must be specified.\n')
}
if(is.null(opt$ld_scores)){
  stop('--ld_scores must be specified.\n')
}
if(is.null(opt$hm3_snplist)){
  stop('--hm3_snplist must be specified.\n')
}
if(is.na(as.numeric(opt$sample_prev))){
  opt$sample_prev<-NULL
}
if(is.na(as.numeric(opt$pop_prev))){
  opt$pop_prev<-NULL
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
log_header(log_file = log_file, opt = opt, script = 'ldsc.R', start.time = start.time)

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

if(opt$hm3_no_mhc & 6 %in% CHROMS){
  # Remove MHC region from hapmap3 SNP-list
  hm3<-fread(opt$hm3_snplist)

  # Read ref pvar
  pvar<-read_pvar(opt$ref_plink_chr, chr=6)

  # Identify SNPs in MHA/HLA region
  # Assumes BP column is GRCh37
  mhc <- pvar[(pvar$CHR == 6 & pvar$BP > 28e6 & pvar$BP < 34e6),]

  # Remove MHC SNPs from hm3 snplist
  hm3_no_mhc<-hm3[!(hm3$SNP %in% mhc$SNP),]

  write.table(hm3_no_mhc, paste0(tmp_dir,'/hm3_no_mhc.snplist'), col.names=T, row.names=F, quote=F)
  opt$hm3_snplist<-paste0(tmp_dir,'/hm3_no_mhc.snplist')

  log_add(log_file = log_file, message = 'MHC region removed for LDSC analysis.')

}

# Read in and write out sumstats to avoid munge error
gwas <- read_sumstats(sumstats = opt$sumstats, chr = CHROMS, log_file = log_file, req_cols = c('CHR','SNP','BP','N','A1','A2','FREQ','BETA','SE','P'))
fwrite(gwas, paste0(tmp_dir,'/sumstats.gz'), row.names=F, quote=F, sep=' ', na='NA')
opt$sumstats<-paste0(tmp_dir,'/sumstats.gz')
GWAS_CHROMS <- unique(gwas$CHR)

ldsc_res <-
  ldsc(
    sumstats = opt$sumstats,
    ldsc = opt$ldsc,
    hm3_snplist = opt$hm3_snplist,
    munge_sumstats = opt$munge_sumstats,
    ld_scores = opt$ld_scores,
    log_file = log_file
  )

write.table(ldsc_res$h2, paste0(opt$output,'.hsq_obs'), col.names=F, row.names = F, quote = F)

# Convert the SNP-heritability to the liability scale
if(!is.null(opt$pop_prev) && !is.null(opt$sample_prev)) {
  ldsc_res$h2_liab <- h2l_R2(k = opt$pop_prev, r2 = ldsc_res$h2, p = opt$sample_prev)
  ldsc_res$h2_liab_se <- h2l_R2(k = opt$pop_prev, r2 = ldsc_res$h2_se, p = opt$sample_prev)
  log_add(log_file = log_file, message = paste0('SNP-heritability estimate on the liability scale = ', round(ldsc_res$h2_liab, 4), " (", round(ldsc_res$h2_liab_se, 4), ")."))
  
  write.table(ldsc_res$h2_liab, paste0(opt$output,'.hsq_liab'), col.names=F, row.names = F, quote = F)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()

