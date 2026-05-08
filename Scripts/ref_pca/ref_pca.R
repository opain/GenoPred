#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--config", action="store", default=NULL, type='character',
      help="Pipeline configuration file [required]"),
  make_option("--population", action="store", default=0.05, type='numeric',
      help="Population to generate PCs [optional]"),
  make_option("--maf", action="store", default=0.05, type='numeric',
      help="Minor allele frequency threshold [optional]"),
  make_option("--geno", action="store", default=0.02, type='numeric',
      help="Variant missingness threshold [optional]"),
  make_option("--hwe", action="store", default=1e-6, type='numeric',
      help="Hardy Weinberg p-value threshold. [optional]"),
  make_option("--n_pcs", action="store", default=6, type='numeric',
      help="Number of PCs [optional]"),
  make_option("--memory", action="store", default=5000, type='numeric',
  		help="Memory limit [optional]"),
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
if(is.null(opt$config)){
  stop('--config must be specified.\n')
}

# Read in outdir
outdir <- read_param(config = opt$config, param = 'outdir', return_obj = F)

# Create output directory
opt$output <- paste0(outdir, '/reference/pc_score_files/', opt$population, '/ref-', opt$population, '-pcs')

# Create output directory
opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste(opt$output,'.log',sep='')
log_header(log_file = log_file, opt = opt, script = 'ref_pca.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

refdir <- read_param(config = opt$config, param = 'refdir', return_obj = F)
if(as.logical(read_param(config = opt$config, param = 'restrict_to_target_variants', return_obj = F))){
  opt$ref_plink_chr <- paste0(outdir, '/reference/ref/ref.chr')
} else {
  opt$ref_plink_chr <- paste0(refdir, '/ref.chr')
}

opt$ref_plink_chr_orig <- opt$ref_plink_chr

if(opt$population != 'TRANS'){
  opt$ref_keep <- paste0(refdir,'/keep_files/', opt$population, '.keep')
} else {
  opt$ref_keep <- NULL
}

opt$plink2 <- 'plink2'

###########
# Extract ref_keep
###########

if(!is.null(opt$ref_keep)){
  plink_subset(keep = opt$ref_keep, chr = CHROMS, plink2 = opt$plink2, pfile = opt$ref_plink_chr, out = paste0(tmp_dir,'/ref_subset.chr'))
  opt$ref_plink_chr<-paste0(tmp_dir,'/ref_subset.chr')
}

###########
# Extract opt$extract
###########

if(!is.null(opt$extract)){
  plink_subset(chr = CHROMS, extract = opt$extract, plink2 = opt$plink2, pfile = opt$ref_plink_chr, out = paste0(tmp_dir,'/ref_subset.chr'))
  opt$ref_plink_chr<-paste0(tmp_dir,'/ref_subset.chr')
}

###########
# QC reference
###########

ref_qc_snplist<-plink_qc_snplist(pfile = opt$ref_plink_chr, chr = CHROMS, plink2 = opt$plink2, geno = opt$geno, maf = opt$maf, hwe = opt$hwe)

###########
# Identify list of LD independent SNPs
###########

log_add(log_file = log_file, message = 'Identifying LD independent SNPs based on reference data.')

# read in reference bim file
ref_pvar<-read_pvar(opt$ref_plink_chr, chr = CHROMS)

# Subset ref_pvar to contain QC'd variants
ref_pvar<-ref_pvar[ref_pvar$SNP %in% ref_qc_snplist,]

# Remove regions of high LD
ref_pvar <- remove_regions(dat = ref_pvar, regions = long_ld_coord)
log_add(log_file = log_file, message = paste0(nrow(ref_pvar),' variants after removal of LD high regions.'))

# Perform LD pruning
ld_indep <- plink_prune(pfile = opt$ref_plink_chr, chr = CHROMS, plink2 = opt$plink2, extract = ref_pvar$SNP)
log_add(log_file = log_file, message = paste0(length(ld_indep),' independent variants retained.'))

###########
# Perform PCA based on reference
###########

log_add(log_file = log_file, message = 'Performing PCA based on reference.')

snp_weights<-plink_pca(pfile = opt$ref_plink_chr, chr = CHROMS, plink2 = opt$plink2, extract = ld_indep, n_pc = opt$n_pcs)
fwrite(snp_weights, paste0(opt$output,'.eigenvec.var'), row.names = F, quote=F, sep=' ', na='NA')

if(file.exists(paste0(opt$output,'.eigenvec.var.gz'))){
  system(paste0('rm ',opt$output,'.eigenvec.var.gz'))
}

system(paste0('gzip ',opt$output,'.eigenvec.var'))

###
# Calculate PCs in the reference sample for scaling the target sample factor scores.
###

log_add(log_file = log_file, message = 'Computing reference PCs.')

# Calculate PCs in the full reference
ref_pcs<-plink_score(pfile = opt$ref_plink_chr_orig, chr = CHROMS, plink2 = opt$plink2, score = paste0(opt$output,'.eigenvec.var.gz'), center = T)

# Save raw PCs
fwrite(ref_pcs, paste0(opt$output, '.profiles'), row.names = F, quote=F, sep=' ', na='NA')

if(opt$population == 'TRANS'){
  ref_pcs_scale <- score_mean_sd(scores = ref_pcs)
  fwrite(ref_pcs_scale, paste0(opt$output, '.TRANS.scale'), row.names = F, quote=F, sep=' ', na='NA')
  scores_scaled<-score_scale(score=ref_pcs, ref_scale=ref_pcs_scale)
  fwrite(scores_scaled, paste0(opt$output, '.profiles'), row.names = F, quote=F, sep=' ', na='NA')
} else {
  # Read reference population data
  pop_data <- read_pop_data(paste0(refdir, '/ref.pop.txt'))
  keep <- pop_data[pop_data$POP == opt$population,]
  ref_pcs_pop <- ref_pcs[paste0(ref_pcs$FID, '_', ref_pcs$IID) %in% paste0(keep$FID, '_', keep$IID),]
  ref_pcs_scale <- score_mean_sd(scores = ref_pcs_pop)
  fwrite(ref_pcs_scale, paste0(opt$output, '.', opt$population, '.scale'), row.names = F, quote=F, sep=' ', na='NA')
  scores_scaled<-score_scale(score=ref_pcs_pop, ref_scale=ref_pcs_scale)
  fwrite(scores_scaled, paste0(opt$output, '.profiles'), row.names = F, quote=F, sep=' ', na='NA')
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
