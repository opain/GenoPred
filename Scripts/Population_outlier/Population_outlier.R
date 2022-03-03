#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--target_plink_chr", action="store", default=NA, type='character',
    help="Path to per chromosome target PLINK files [required]"),
make_option("--target_plink", action="store", default=NA, type='character',
    help="Path to per chromosome target PLINK files [required]"),
make_option("--target_fam", action="store", default=NA, type='character',
    help="Target sample fam file. [optional]"),
make_option("--maf", action="store", default=NA, type='numeric',
    help="Minor allele frequency threshold [optional]"),
make_option("--geno", action="store", default=0.02, type='numeric',
    help="Variant missingness threshold [optional]"),
make_option("--hwe", action="store", default=NA, type='numeric',
    help="Hardy Weinberg p-value threshold. [optional]"),
make_option("--n_pcs", action="store", default=10, type='numeric',
		help="Number of PCs (min=4) [optional]"),
make_option("--plink", action="store", default='plink', type='character',
		help="Path PLINK software binary [required]"),
make_option("--plink2", action="store", default='plink', type='character',
		help="Path PLINK software binary [required]"),
make_option("--target_keep", action="store", default=NA, type='character',
    help="File containing list of keep files [optional]"),
make_option("--output", action="store", default='./PC_projector_output/Output', type='character',
    help="Path for output files [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(ggplot2)
library(cowplot)
library(NbClust)
library(GGally)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# Population_outlier.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

if(!is.na(opt$target_keep)){
  keep_list<-fread(opt$target_keep, header=F)
  keep_list$pop<-gsub('\\..*','',gsub('.*model_pred.','',gsub('.*/','', keep_list$V1)))
  names(keep_list)<-c('file','pop')
} else {
  keep_list<-data.frame(file=NA,
                        pop='')
}

############
# Create file listing variants in regions of long range LD
############

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Excluding long range LD...')
sink()

# Read in the bim file
if(is.na(opt$target_plink)){
  targ_bim<-NULL
  for(i in 1:22){
    targ_bim<-rbind(targ_bim, fread(paste0(opt$target_plink_chr,i,'.bim')))
  }
} else {
  targ_bim<-fread(paste0(opt$target_plink,'.bim'))
}

# Create file removing these regions.
long_ld_exclude<-targ_bim$V2[ (targ_bim$V1 == 1 & targ_bim$V4 >= 48e6 & targ_bim$V4 <= 52e6) |
                               (targ_bim$V1 == 2 & targ_bim$V4 >= 86e6 & targ_bim$V4 <= 100.5e6) |
                               (targ_bim$V1 == 2 & targ_bim$V4 >= 134.5e6 & targ_bim$V4 <= 138e6) |
                               (targ_bim$V1 == 2 & targ_bim$V4 >= 183e6 & targ_bim$V4 <= 190e6) |
                               (targ_bim$V1 == 3 & targ_bim$V4 >= 47.5e6 & targ_bim$V4 <= 50e6) |
                               (targ_bim$V1 == 3 & targ_bim$V4 >= 83.5e6 & targ_bim$V4 <= 87e6) |
                               (targ_bim$V1 == 3 & targ_bim$V4 >= 89e6 & targ_bim$V4 <= 97.5e6) |
                               (targ_bim$V1 == 5 & targ_bim$V4 >= 44.5e6 & targ_bim$V4 <= 50.5e6) |
                               (targ_bim$V1 == 5 & targ_bim$V4 >= 98e6 & targ_bim$V4 <= 100.5e6) |
                               (targ_bim$V1 == 5 & targ_bim$V4 >= 129e6 & targ_bim$V4 <= 132e6) |
                               (targ_bim$V1 == 5 & targ_bim$V4 >= 135.5e6 & targ_bim$V4 <= 138.5e6) |
                               (targ_bim$V1 == 6 & targ_bim$V4 >= 25.5e6 & targ_bim$V4 <= 33.5e6) |
                               (targ_bim$V1 == 6 & targ_bim$V4 >= 57e6 & targ_bim$V4 <= 64e6) |
                               (targ_bim$V1 == 6 & targ_bim$V4 >= 140e6 & targ_bim$V4 <= 142.5e6) |
                               (targ_bim$V1 == 7 & targ_bim$V4 >= 55e6 & targ_bim$V4 <= 66e6) |
                               (targ_bim$V1 == 8 & targ_bim$V4 >= 8e6 & targ_bim$V4 <= 12e6) |
                               (targ_bim$V1 == 8 & targ_bim$V4 >= 43e6 & targ_bim$V4 <= 50e6) |
                               (targ_bim$V1 == 8 & targ_bim$V4 >= 112e6 & targ_bim$V4 <= 115e6) |
                               (targ_bim$V1 == 10 & targ_bim$V4 >= 37e6 & targ_bim$V4 <= 43e6) |
                               (targ_bim$V1 == 11 & targ_bim$V4 >= 46e6 & targ_bim$V4 <= 57e6) |
                               (targ_bim$V1 == 11 & targ_bim$V4 >= 87.5e6 & targ_bim$V4 <= 90.5e6) |
                               (targ_bim$V1 == 12 & targ_bim$V4 >= 33e6 & targ_bim$V4 <= 40e6) |
                               (targ_bim$V1 == 12 & targ_bim$V4 >= 109.5e6 & targ_bim$V4 <= 112e6) |
                               (targ_bim$V1 == 20 & targ_bim$V4 >= 32e6 & targ_bim$V4 <= 34.5e6)]

write.table(long_ld_exclude, paste0(opt$output_dir,'long_ld.exclude'), col.names=F, row.names=F, quote=F)

rm(long_ld_exclude,targ_bim)
gc()

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

########
# Within ancestry QC
########
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('QC qill be performed for', paste(keep_list$pop, collapse=', '),'populations...\n')
sink()

for(pop in keep_list$pop){

  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('~~~~~~  ', pop,'    ~~~~~~\n')
  sink()
  
  ###########
  # Perform PCA on QC'd and independent variants
  ###########

  # Identify LD independent QC'd SNPs using a random subset 1000 individuals
  keep_file<-fread(keep_list$file[keep_list$pop == pop],header=FALSE)
  if(nrow(keep_file) < 50){
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Skipped due to insufficient sample size.\n')
    cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    sink()
      
    write.table(keep_file, paste0(opt$output,'.',pop,'.keep'), col.names=F, row.names=F, quote=F)
    rm(keep_file)
    next
  }

  if (nrow(keep_file) < 100){
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat(paste0('Warning: sample size for ancestry ',pop,' is smaller than 100 . (',nrow(keep_file),') \n'))
    sink()
  }
  
  keep_file_subset<-keep_file[sample(min(1000, nrow(keep_file))),]
  write.table(keep_file_subset, paste0(opt$output_dir,pop,'_subset.keep'), col.names=F, row.names=F, quote=F)
  
  rm(keep_file,keep_file_subset)
  gc()
  
  if(is.na(opt$target_fam)){
    if(is.na(opt$target_plink)){
      for(i in 1:22){
        system(paste0(opt$plink,' --bfile ',opt$target_plink_chr,i,' --keep ',opt$output_dir,pop,'_subset.keep --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --exclude ',opt$output_dir,'long_ld.exclude --indep-pairwise 1000 5 0.2 --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
      }
    } else {
      system(paste0(opt$plink,' --bfile ',opt$target_plink,' --keep ',opt$output_dir,pop,'_subset.keep --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --exclude ',opt$output_dir,'long_ld.exclude --indep-pairwise 1000 5 0.2 --chr 1-22 --write-snplist --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
    }
  } else {
    if(is.na(opt$target_plink)){
      for(i in 1:22){
        system(paste0(opt$plink,' --bfile ',opt$target_plink_chr,i,' --fam ',opt$target_fam,' --keep ',opt$output_dir,pop,'_subset.keep --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --exclude ',opt$output_dir,'long_ld.exclude --indep-pairwise 1000 5 0.2 --write-snplist --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
      }
    } else {
      system(paste0(opt$plink,' --bfile ',opt$target_plink,' --fam ',opt$target_fam,' --keep ',opt$output_dir,pop,'_subset.keep --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --exclude ',opt$output_dir,'long_ld.exclude --indep-pairwise 1000 5 0.2 --chr 1-22 --write-snplist --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
    }
  }

  # Randomly prune to 100,000 variants
  if(is.na(opt$target_plink)){
    prune_all<-NULL
    for(i in 1:22){
      prune_all<-rbind(prune_all, fread(paste0(opt$output_dir,'target.QC.chr',i,'.prune.in'), header=F))
    }
    write.table(prune_all, paste0(opt$output_dir,'target.QC.prune.in'), col.names=F, row.names=F, quote=F)
    
  }
  
  # Calculate SNP weights
  if(is.na(opt$target_fam)){
    if(is.na(opt$target_plink)){
      # Create merge list
      ref_merge_list<-paste0(opt$target_plink_chr,1:22)
      write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)
      
      system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --keep ',opt$output_dir,pop,'_subset.keep --threads 1  --pca ',opt$n_pcs,' var-wts --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
    } else {
      system(paste0(opt$plink,' --bfile ',opt$target_plink,' --keep ',opt$output_dir,pop,'_subset.keep --threads 1 --pca ',opt$n_pcs,' var-wts --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
    }
  } else {
    if(is.na(opt$target_plink)){
      # Create merge list
      ref_merge_list<-paste0(opt$target_plink_chr,1:22)
      write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)
      
      system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --fam ',opt$target_fam,' --keep ',opt$output_dir,pop,'_subset.keep --threads 1  --pca ',opt$n_pcs,' var-wts --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
    } else {
      system(paste0(opt$plink,' --bfile ',opt$target_plink,' --fam ',opt$target_fam,' --keep ',opt$output_dir,pop,'_subset.keep --threads 1 --pca ',opt$n_pcs,' var-wts --extract ',opt$output_dir,'target.QC.prune.in --chr 1-22 --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
    }
  }
  
  # Modify .eigenvec.var file to contain header (not sure why there is no header)
  var_file<-fread(paste0(opt$output_dir,'target.QC.eigenvec.var'), header=F)
  names(var_file)<-c('CHR','SNP','A1','A2',paste0('PC',1:opt$n_pcs))
  write.table(var_file, paste0(opt$output_dir,'target.QC.eigenvec.var'), col.names=T, row.names=F, quote=F)
  
  rm(var_file)
  gc()

  # Calculate PCs in the full sample
  system(paste0("cut -f 2 -d ' ' ",opt$output_dir,'target.QC.eigenvec.var | tail -n +2 > ',opt$output_dir,'score_file.snplist'))
  
  if(is.na(opt$target_fam)){
    if(is.na(opt$target_plink)){
      for(i in 1:22){
        system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --extract ',opt$output_dir,'score_file.snplist --keep ',keep_list$file[keep_list$pop == pop],' --score ',opt$output_dir,'target.QC.eigenvec.var header-read 2 3 no-mean-imputation --threads 1 --score-col-nums 5-',as.numeric(opt$n_pcs)+4,' --out ',opt$output_dir,'profiles.chr',i,' --memory ',floor(opt$memory*0.9)))
      }
    } else {
      system(paste0(opt$plink2, ' --bfile ',opt$target_plink,' --extract ',opt$output_dir,'score_file.snplist --keep ',keep_list$file[keep_list$pop == pop],' --score ',opt$output_dir,'target.QC.eigenvec.var header-read 2 3 no-mean-imputation --threads 1 --score-col-nums 5-',as.numeric(opt$n_pcs)+4,' --out ',opt$output_dir,'profiles',i,' --memory ',floor(opt$memory*0.9)))
    }
  } else {
    if(is.na(opt$target_plink)){
      for(i in 1:22){
        system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --fam ',opt$target_fam,' --extract ',opt$output_dir,'score_file.snplist --keep ',keep_list$file[keep_list$pop == pop],' --score ',opt$output_dir,'target.QC.eigenvec.var header-read 2 3 no-mean-imputation --threads 1 --score-col-nums 5-',as.numeric(opt$n_pcs)+4,' --out ',opt$output_dir,'profiles.chr',i,' --memory ',floor(opt$memory*0.9)))
      }
    } else {
      system(paste0(opt$plink2, ' --bfile ',opt$target_plink,' --fam ',opt$target_fam,' --extract ',opt$output_dir,'score_file.snplist --keep ',keep_list$file[keep_list$pop == pop],' --score ',opt$output_dir,'target.QC.eigenvec.var header-read 2 3 no-mean-imputation --threads 1 --score-col-nums 5-',as.numeric(opt$n_pcs)+4,' --out ',opt$output_dir,'profiles --memory ',floor(opt$memory*0.9)))
    }
  }
  
  if(is.na(opt$target_plink)){
    # Add up the scores across chromosomes
    scores<-fread(cmd=paste0('cut -f 1-2 ',opt$output_dir,'profiles.chr22.sscore'))
    names(scores)<-c('FID','IID')
    
    var_list<-fread(paste0(opt$output_dir,'target.QC.eigenvec.var'))
    nsnp_all<-0
    for(i in 1:22){
    	sink(file = paste(opt$output,'.log',sep=''), append = T)
    	cat('Adding up PC projections for chromosome ', i, ', ancestry ',pop,' ... \n')
    	sink()
    	profile<-data.frame(fread(paste0(opt$output_dir,'profiles.chr',i,'.sscore')))
    	profile<-as.matrix(profile[,grepl('PC',names(profile))])
    	bim<-fread(paste0(opt$target_plink_chr,i,'.bim'))
    	nsnp<-sum(bim$V2 %in% var_list$SNP)
    	nsnp_all<-nsnp_all+nsnp
    	profile<-profile*nsnp
    	if(i == 1){
    		profile_all<-profile
    	} else {
    		profile_all<-profile_all+profile	
    	}
    }
    
    profile_all<-profile_all/nsnp_all
  } else {
    scores<-fread(cmd=paste0('cut -f 1-2 ',opt$output_dir,'profiles.sscore'))
    names(scores)<-c('FID','IID')
  
    profile_all<-data.frame(fread(paste0(opt$output_dir,'profiles.sscore')))
  
    profile_all<-as.matrix(profile_all[,grepl('PC',names(profile_all))])
  }
  
  profile_all<-data.table(profile_all)
  names(profile_all)<-paste0('PC',1:as.numeric(opt$n_pcs))
  
  scores<-cbind(scores, profile_all)
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Done!\n')
  sink()
  
  targ_PCs<-data.frame(scores)
  targ_PCs[,-1:-2]<-scale(targ_PCs[,-1:-2])
  rm(scores,profile_all,var_list,nsnp_all)
  gc()
  
  ###
  # Clean up temporary files
  ###
  
  system(paste0('rm ',opt$output_dir,'profiles*'))

  ###
  # Create plot PC scores of target sample
  ###
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Plotting target sample PCs...\n')
  sink()
  
  pairs_plot<-ggpairs(targ_PCs[,-1:-2])

  png(paste0(opt$output,'.',pop,'.PCs_plot.png'), units='px', res=300, width=6000, height=6000)
  print(pairs_plot)
  dev.off()
  
  rm(pairs_plot)
  gc()
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Done!\n')
  sink()
  
  ###
  # Identify outliers
  ###
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Defining number of centroids...\n')
  sink()
  
  png(paste0(opt$output,'.',pop,'.NbClust.png'), units='px', res=300, width=3000, height=3000)
  # Define centroids
  if(dim(targ_PCs)[1] > 1000){
  	n_clust_sol<-NbClust(data = targ_PCs[sample(1000),-1:-2], distance = "euclidean", min.nc = 2, max.nc = 10, method = 'kmeans', index='all')
  } else {
  	n_clust_sol<-NbClust(data = targ_PCs[,-1:-2], distance = "euclidean", min.nc = 2, max.nc = 10, method = 'kmeans', index='all')
  }
  dev.off()
  
  n_clust_opt<-length(unique(n_clust_sol$Best.partition))
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat(n_clust_opt, 'centroids was found to fit best.\n')
  sink()
  
  rm(n_clust_sol)
  gc()

  # Calculate distance from centroids
  k_res<-kmeans(targ_PCs[-1:-2], n_clust_opt)
  k_res$centers<-data.frame(k_res$centers)
  
  targ_PCs_distance<-list()
  pc_var <- sapply(targ_PCs[,3:(opt$n_pcs+2)], var)
    
  for(centroid_n in 1:n_clust_opt){
    targ_PCs_distance[[centroid_n]]<-data.frame(targ_PCs[i,])
    targ_PCs_distance[[centroid_n]][,3:(opt$n_pcs + 2)] <- t(t(as.matrix(targ_PCs_distance[[centroid_n]][3:(opt$n_pc+2)])) - unlist(k_res$centers[centroid_n,]))^2 # calculate distance
    targ_PCs_distance[[centroid_n]][,3:(opt$n_pcs + 2)] <- targ_PCs_distance[[centroid_n]][,3:(opt$n_pcs + 2)] / pc_var # scale by variance (variance is 1. because we scaled the PCs before...)
  }
  
  rm(k_res, pc_var)
  gc()

  # Plot centroid distance
  pdf(paste0(opt$output,'.',pop,'.centroid_distance_plot.pdf'))
  for(centroid_n in 1:n_clust_opt){
    for(pc_n in 1:opt$n_pcs){
      hist(targ_PCs_distance[[centroid_n]][[paste0('PC',pc_n)]], main=paste0('Centroid ', centroid_n,', PC', pc_n))
    }
  }
  dev.off()
  
  # Remove outliers using a cut off defined as Q3+30*IQR
  outliers<-NULL
  for(centroid_n in 1:n_clust_opt){
    for(pc_n in 1:opt$n_pcs){
      
      Q <- quantile(targ_PCs_distance[[centroid_n]][[paste0('PC',pc_n)]], probs=c(.75), na.rm = FALSE)
      IQR <- IQR(targ_PCs_distance[[centroid_n]][[paste0('PC',pc_n)]])
      
      outliers<-c(outliers,which(targ_PCs_distance[[centroid_n]][[paste0('PC',pc_n)]] > Q+30*IQR))
    }
  }
  outliers<-unique(outliers)
    
  if (length(outliers) > 0){
      
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Found ',length(outliers), ' for population ',pop,'.\n')
    sink()
      
    targ_PCs_distance_no_outliers<-targ_PCs_distance
    for(centroid_n in 1:n_clust_opt){
      for(pc_n in 1:opt$n_pcs){
        targ_PCs_distance_no_outliers[[centroid_n]][[paste0('PC',pc_n)]][outliers]<-NA
      }
      targ_PCs_distance_no_outliers[[centroid_n]]<-targ_PCs_distance_no_outliers[[centroid_n]][complete.cases(targ_PCs_distance_no_outliers[[centroid_n]]),]
    }
    
    # Plot centroid distance afte removing outliers
    pdf(paste0(opt$output,'.',pop,'.centroid_distance_noOutlier_plot.pdf'))
    for(centroid_n in 1:n_clust_opt){
      for(pc_n in 1:opt$n_pcs){
        hist(targ_PCs_distance_no_outliers[[centroid_n]][[paste0('PC',pc_n)]], main=paste0('Centroid ', centroid_n,', PC', pc_n))
      }
    }
    dev.off()
    
    rm(targ_PCs_distance_no_outliers)
    gc()
      
    targ_PCs_noOutliers<-targ_PCs[-outliers,]
      
    pairs_plot<-ggpairs(targ_PCs_noOutliers[,-1:-2])

    png(paste0(opt$output,'.',pop,'.PCs_noOutlier_plot.png'), units='px', res=300, width=6000, height=6000)
    print(pairs_plot)
    dev.off()
  
    rm(pairs_plot)
    gc()
    
  } else {
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Found no outliers for population ',pop,'.\n')
    sink()
      
    targ_PCs_noOutliers<-targ_PCs
  }
      
  ###
  # Create plot PC scores of target sample after removal of outliers
  ###
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat(dim(targ_PCs_noOutliers)[1],'out of', dim(targ_PCs)[1],'remain.\n')
  sink()
  
  ###
  # Save PCs, distances and keep file for each population
  ###
  
  write.table(targ_PCs_noOutliers[,c('FID','IID')], paste0(opt$output,'.',pop,'.keep'), col.names=F, row.names=F, quote=F)
  fwrite(targ_PCs, paste0(opt$output,'.',pop,'.PCs.txt.gz'), col.names=F, row.names=F, quote=F)
    
  saveRDS(targ_PCs_distance, file=paste0(opt$output,'.',pop,'.centroid_distance.RDS'))

  ###
  # clean up unwanted files
  ###
  
  rm(targ_PCs_noOutliers,targ_PCs,outliers,targ_PCs_distance)
  gc()

  system(paste0('rm ',opt$output_dir, 'target*'))
  system(paste0('rm ',opt$output_dir, 'score_file.snplist*'))
  system(paste0('rm ',opt$output_dir, pop,'_subset.keep'))
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
  sink()
  
}

system(paste0('rm ',opt$output_dir, 'long_ld.exclude'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
