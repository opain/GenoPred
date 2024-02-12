#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--ref_plink_chr", action="store", default=NULL, type='character',
		help="Path to per chromosome reference PLINK2 files [required]"),
make_option("--ref_keep", action="store", default=NULL, type='character',
		help="Keep file to subset individuals in reference for PCA [required]"),
make_option("--maf", action="store", default=0.05, type='numeric',
    help="Minor allele frequency threshold [optional]"),
make_option("--geno", action="store", default=0.02, type='numeric',
    help="Variant missingness threshold [optional]"),
make_option("--hwe", action="store", default=1e-6, type='numeric',
    help="Hardy Weinberg p-value threshold. [optional]"),
make_option("--n_pcs", action="store", default=6, type='numeric',
		help="Number of PCs [optional]"),
make_option("--plink2", action="store", default='plink2', type='character',
		help="Path PLINKv2 software binary [optional]"),
make_option("--output", action="store", default=NULL, type='character',
		help="Path for output files [required]"),
make_option("--pop_data", action="store", default=NULL, type='character',
    help="Population data for the reference samples [required]"),
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
if(is.null(opt$pop_data)){
  stop('--pop_data must be specified.\n')
}
if(is.null(opt$ref_plink_chr)){
  stop('--ref_plink_chr must be specified.\n')
}
if(is.null(opt$output)){
  stop('--output must be specified.\n')
}

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

###########
# Extract ref_keep
###########

if(!is.null(opt$ref_keep)){
  plink_subset(keep = opt$ref_keep, chr = CHROMS, plink2 = opt$plink2, pfile = opt$ref_plink_chr, out = paste0(tmp_dir,'/ref_subset.chr'))
  opt$ref_plink_chr_subset<-paste0(tmp_dir,'/ref_subset.chr')
} else {
  opt$ref_plink_chr_subset<-opt$ref_plink_chr
}

###########
# QC reference
###########

ref_qc_snplist<-plink_qc_snplist(pfile = opt$ref_plink_chr_subset, chr = CHROMS, plink2 = opt$plink2, geno = opt$geno, maf = opt$maf, hwe = opt$hwe)

###########
# Identify list of LD independent SNPs
###########

log_add(log_file = log_file, message = 'Identifying LD independent SNPs based on reference data.')

# read in reference bim file
ref_pvar<-read_pvar(opt$ref_plink_chr_subset, chr = CHROMS)

# Subset ref_bim to contain QC'd variants
ref_pvar<-ref_pvar[ref_pvar$SNP %in% ref_qc_snplist,]

# Remove regions of high LD
ref_pvar <- remove_regions(dat = ref_pvar, regions = long_ld_coord)
log_add(log_file = log_file, message = paste0(nrow(ref_pvar),' variants after removal of LD high regions.'))

# Perform LD pruning
ld_indep <- plink_prune(pfile = opt$ref_plink_chr_subset, chr = CHROMS, plink2 = opt$plink2, extract = ref_bim$SNP)
log_add(log_file = log_file, message = paste0(length(ld_indep),' independent variants retained.'))

###########
# Perform PCA based on reference
###########

log_add(log_file = log_file, message = 'Performing PCA based on reference.')

snp_weights<-plink_pca(pfile = opt$ref_plink_chr_subset, chr = CHROMS, plink2 = opt$plink2, extract = ld_indep, n_pc = opt$n_pcs)
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
ref_pcs<-plink_score(pfile = opt$ref_plink_chr, chr = CHROMS, plink2 = opt$plink2, score = paste0(opt$output,'.eigenvec.var.gz'))

# Calculate scale within each reference population
pop_data<-fread(opt$pop_data)
pop_data<-data.table(
  FID=pop_data$`#IID`,
  IID=pop_data$`#IID`,
  POP=pop_data$POP
)

for(pop_i in unique(pop_data$POP)){
  ref_pcs_scale_i <- score_mean_sd(scores = ref_pcs, keep = pop_data[pop_data$POP == pop_i, c('FID','IID'), with=F])
  fwrite(ref_pcs_scale_i, paste0(opt$output, '.', pop_i, '.scale'), row.names = F, quote=F, sep=' ', na='NA')
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
