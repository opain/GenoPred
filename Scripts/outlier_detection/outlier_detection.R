#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
  make_option("--target_plink_chr", action="store", default=NULL, type='character',
              help="Path to per chromosome target PLINK files [required]"),
  make_option("--maf", action="store", default=0.05, type='numeric',
              help="Minor allele frequency threshold [optional]"),
  make_option("--geno", action="store", default=0.02, type='numeric',
              help="Variant missingness threshold [optional]"),
  make_option("--hwe", action="store", default=1e-6, type='numeric',
              help="Hardy Weinberg p-value threshold. [optional]"),
  make_option("--n_pcs", action="store", default=10, type='numeric',
              help="Number of PCs (min=4) [optional]"),
  make_option("--plink2", action="store", default='plink2', type='character',
              help="Path PLINK2 software binary [required]"),
  make_option("--keep_list", action="store", default=NULL, type='character',
              help="File containing list of keep files and corresponding population code [optional]"),
  make_option("--unrel", action="store", default=NA, type='character',
              help="File containing list of unrelated individuals [optional]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
              help="Number of cores for parallel computing [optional]"),
  make_option("--test", action="store", default=NA, type='character',
              help="Specify test mode [optional]"),
  make_option("--output", action="store", default=NULL, type='character',
              help="Path for output files [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')
library(ggplot2)
library(cowplot)
library(NbClust)
library(GGally)

set.seed(1)

# Check required inputs
if(is.null(opt$target_plink_chr)){
  stop('--target_plink_chr must be specified.\n')
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
log_header(log_file = log_file, opt = opt, script = 'outlier_detection.R', start.time = start.time)

if(!is.na(opt$unrel) && opt$unrel == 'NA'){
  opt$unrel<-NA
}

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

######
# Estimate relatedness
######

if(is.na(opt$unrel)){
  log_add(log_file = log_file, message = 'Estimating relatedness.')

  # Identify high quality variants
  target_qc_snplist<-plink_qc_snplist(pfile = opt$target_plink_chr, chr = CHROMS, plink2 = opt$plink2, geno = opt$geno, maf = opt$maf, hwe = opt$hwe)

  # Generate kinship matrix and list of unrelated individuals
  plink_king(pfile = opt$target_plink_chr, chr = CHROMS, extract = target_qc_snplist, plink2 = opt$plink2, out = opt$output, threads = opt$n_cores)

  # Read in list of unrelated individuals
  unrel <- fread(paste0(opt$output, '.unrelated.keep'), header=F)
  names(unrel)<-c('FID','IID')
} else {
  log_add(log_file = log_file, message = 'Using user-specified list of unrelated individuals.')

  # Read in list of unrelated individuals
  unrel <- fread(opt$unrel, header=F)
  if(ncol(unrel) == 1){
    if(names(unrel)[1] == 'V1'){
      unrel$V2<-unrel$V1
    }
  }
  unrel<-unrel[, 1:2, with=F]
  names(unrel)[1:2]<-c('FID','IID')
}

log_add(log_file = log_file, message = paste0('There are ', nrow(unrel), ' unrelated individuals present.'))

###
# Read in keep_list
###

if(!is.null(opt$keep_list)){
  keep_list <- fread(opt$keep_list)
} else {
  keep_list<-data.frame(POP='',
                        file=NA)
}

############
# Create file listing variants in regions of long range LD
############

log_add(log_file = log_file, message = 'Excluding long range LD regions.')

targ_pvar <- read_pvar(opt$target_plink_chr, chr = CHROMS)
targ_pvar <- remove_regions(dat = targ_pvar, regions = long_ld_coord)

log_add(log_file = log_file, message = paste0(nrow(targ_pvar),' variants after removal of LD high regions.'))

########
# Within ancestry QC
########

if(!is.null(opt$keep_list)){
  log_add(log_file = log_file, message = paste0('QC will be performed for ', paste(keep_list$POP, collapse=', '),' populations.'))
}

for(pop in keep_list$POP){

  log_add(log_file = log_file, message = paste0('~~~~~~  ', pop,'    ~~~~~~'))

  # Read in keep file for population
  keep_file <- fread(keep_list$file[keep_list$POP == pop], header=F)
  if(ncol(keep_file) == 1){
    keep_file <- data.table(
      FID = keep_file$V1,
      IID = keep_file$V1)
  } else {
    keep_file <- data.table(
      FID = keep_file$V1,
      IID = keep_file$V2)
  }

  if(nrow(keep_file) < 100){
    log_add(log_file = log_file, message = c('Skipped due to insufficient sample size.', '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'))
    next
  }

  ###########
  # Perform PCA on QC'd and independent variants
  ###########

  # Create QC'd SNP-list
  target_qc_snplist <- plink_qc_snplist(pfile = opt$target_plink_chr, plink2 = opt$plink2, chr = CHROMS, keep = keep_file, maf = opt$maf, geno = opt$geno, hwe = opt$hwe, threads = opt$n_cores)

  # Remove high LD regions
  target_qc_snplist <- target_qc_snplist[target_qc_snplist %in% targ_pvar$SNP]

  # Perform LD pruning
  ld_indep <- plink_prune(pfile = opt$target_plink_chr, chr = CHROMS, keep = keep_file, plink2 = opt$plink2, extract = target_qc_snplist, threads = opt$n_cores)

  # To improve efficiency, derive PCs using random subset of 1000 individuals.
  # These individuals should be unrelated
  keep_file_subset <- merge(keep_file, unrel, by=c('FID','IID'))
  keep_file_subset <- keep_file_subset[sample(1000, replace = T),]
  keep_file_subset <- keep_file_subset[!duplicated(keep_file_subset),]

  # Run PCA
  snp_weights <- plink_pca(pfile = opt$target_plink_chr, keep = keep_file_subset, chr = CHROMS, plink2 = opt$plink2, extract = ld_indep, n_pc = opt$n_pcs, threads = opt$n_cores)
  fwrite(snp_weights, paste0(tmp_dir,'/ref.eigenvec.var'), row.names = F, quote=F, sep=' ', na='NA')

  # Project into the full population
  target_pcs <- plink_score(pfile = opt$target_plink_chr, keep = keep_file, chr = CHROMS, plink2 = opt$plink2, score = paste0(tmp_dir,'/ref.eigenvec.var'), threads = opt$n_cores)

  # Create plot PC scores of target sample
  pairs_plot <- ggpairs(target_pcs[,-1:-2])

  png(paste0(opt$output, '.', pop, '.pc_plot.png'), units = 'px', res = 300, width = 6000, height = 6000)
  print(pairs_plot)
  dev.off()

  ###
  # Identify outliers
  ###

  log_add(log_file = log_file, message = 'Defining number of centroids')

  png(paste0(opt$output, '.', pop, '.NbClust.png'), units = 'px', res = 300, width = 3000, height = 2000)
  # Define number of centroids, again using a random subset of 1000 unrelated individuals
  target_pcs_subset <- merge(target_pcs, keep_file_subset, by.x = c('FID','IID'), by.y = c('FID','IID'))
  n_clust_sol <- NbClust(data = target_pcs_subset[,-1:-2], distance = "euclidean", min.nc = 2, max.nc = 10, method = 'kmeans', index='all')
  dev.off()

  n_clust_opt<-length(unique(n_clust_sol$Best.partition))

  log_add(log_file = log_file, message = paste0(n_clust_opt, ' centroids was found to fit best.'))

  # Identify centroids using kmeans and calculate distance from them
  k_res <- kmeans(target_pcs[,-1:-2], n_clust_opt)
  k_res$centers<-data.frame(k_res$centers)

  target_pcs_distance<-list()
  for(centroid_n in 1:n_clust_opt){
    target_pcs_distance[[centroid_n]] <- data.frame(target_pcs[, c('FID', 'IID'), with = F])
    for(pc_n in 1:opt$n_pcs){
      # Subtract the mean of centroid from PC
      target_pcs_distance[[centroid_n]][paste0('PC',pc_n)]<-k_res$centers[[paste0('PC',pc_n)]][centroid_n]-target_pcs[[paste0('PC',pc_n)]]

      # square the value (to make all positive)
      target_pcs_distance[[centroid_n]][paste0('PC',pc_n)]<-target_pcs_distance[[centroid_n]][paste0('PC',pc_n)]^2

      # Divide by the variance of the PC
      target_pcs_distance[[centroid_n]][paste0('PC',pc_n)]<-target_pcs_distance[[centroid_n]][paste0('PC',pc_n)]/var(target_pcs[[paste0('PC',pc_n)]])

      # Scale the distance from centroids
      target_pcs_distance[[centroid_n]][paste0('PC',pc_n)]<-target_pcs_distance[[centroid_n]][paste0('PC',pc_n)]/var(target_pcs[[paste0('PC',pc_n)]])

    }
  }

  # Plot centroid distance
  pdf(paste0(opt$output, '.', pop, '.centroid_distance_plot.pdf'))
  for(centroid_n in 1:n_clust_opt){
    for(pc_n in 1:opt$n_pcs){
      hist(target_pcs_distance[[centroid_n]][[paste0('PC', pc_n)]], main = paste0('Centroid ', centroid_n,', PC', pc_n))
    }
  }
  dev.off()

  # Remove outliers using a cut off defined as Q3+30*IQR
  outliers<-NULL
  for(centroid_n in 1:n_clust_opt){
    for(pc_n in 1:opt$n_pcs){

      Q <- quantile(target_pcs_distance[[centroid_n]][[paste0('PC',pc_n)]], probs=c(.75), na.rm = FALSE)
      IQR <- IQR(target_pcs_distance[[centroid_n]][[paste0('PC',pc_n)]])

      outliers<-c(outliers,which(target_pcs_distance[[centroid_n]][[paste0('PC',pc_n)]] > Q+30*IQR))
    }
  }
  outliers<-unique(outliers)

  target_pcs_distance_no_outliers<-target_pcs_distance
  for(centroid_n in 1:n_clust_opt){
    for(pc_n in 1:opt$n_pcs){
      target_pcs_distance_no_outliers[[centroid_n]][[paste0('PC',pc_n)]][outliers]<-NA
    }
    target_pcs_distance_no_outliers[[centroid_n]]<-target_pcs_distance_no_outliers[[centroid_n]][complete.cases(target_pcs_distance_no_outliers[[centroid_n]]),]
  }

  # Plot centroid distance afte removing outliers
  pdf(paste0(opt$output, '.', pop, '.centroid_distance_noOutlier_plot.pdf'))
  for(centroid_n in 1:n_clust_opt){
    for(pc_n in 1:opt$n_pcs){
      hist(target_pcs_distance_no_outliers[[centroid_n]][[paste0('PC',pc_n)]], main = paste0('Centroid ', centroid_n,', PC', pc_n))
    }
  }
  dev.off()

  ###
  # Create plot PC scores of target sample after removal of outliers
  ###

  target_pcs_noOutliers<-target_pcs[-outliers,]

  pairs_plot<-ggpairs(target_pcs_noOutliers[,-1:-2])

  png(paste0(opt$output,'.',pop,'.pc_noOutlier_plot.png'), units='px', res=300, width=6000, height=6000)
  print(pairs_plot)
  dev.off()

  log_add(log_file = log_file, message = paste0(nrow(target_pcs_noOutliers),' out of ', nrow(target_pcs),' remain.'))

  ###
  # Save PCs, distances and keep file for each population
  ###

  write.table(target_pcs_noOutliers[,c('FID', 'IID')], paste0(opt$output, '.', pop, '.keep'), col.names = F, row.names = F, quote = F)
  write.table(target_pcs, paste0(opt$output, '.', pop, '.PCs.txt'), col.names = T, row.names = F, quote = F)
  saveRDS(target_pcs_distance, file = paste0(opt$output, '.', pop, '.centroid_distance.RDS'))

  log_add(log_file = log_file, message = '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()
