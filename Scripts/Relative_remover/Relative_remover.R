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
make_option("--plink", action="store", default='plink', type='character',
		help="Path PLINK software binary [required]"),
make_option("--target_keep", action="store", default=NA, type='character',
    help="File containing list of keep files [optional]"),
make_option("--rel_thresh", action="store", default=0.044, type='numeric',
    help="Kinship relatedness threshold [optional]"),
make_option("--output", action="store", default='./PC_projector_output/Output', type='character',
    help="Path for output files [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(ggplot2)
library(cowplot)

opt$output_dir<-dirname(opt$output)
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# Relative_remover.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

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

###########
# Perform PCA on QC'd and independent variants
###########

# Identify LD independent QC'd SNPs using a random subset 1000 individuals
if(is.na(opt$target_keep)){
  if(is.na(opt$target_fam)){
    if(is.na(opt$target_plink)){
      fam_file<-fread(paste0(opt$target_plink_chr,'22.fam'))
      fam_file<-fam_file[sample(1000),]
      write.table(fam_file, paste0(opt$output_dir,'_subset.keep'), col.names=F, row.names=F, quote=F)
    } else {
      fam_file<-fread(paste0(opt$target_plink,'.fam'))
      fam_file<-fam_file[sample(1000),]
      write.table(fam_file, paste0(opt$output_dir,'_subset.keep'), col.names=F, row.names=F, quote=F)
    }
      fam_file<-fread(paste0(opt$target_plink,'.fam'))
      fam_file<-fam_file[sample(1000),]
      write.table(fam_file, paste0(opt$output_dir,'_subset.keep'), col.names=F, row.names=F, quote=F)
  } else {
      fam_file<-fread(opt$target_fam)
      fam_file<-fam_file[sample(1000),]
      write.table(fam_file, paste0(opt$output_dir,'_subset.keep'), col.names=F, row.names=F, quote=F)
  } 

} else {

  keep_file<-fread(opt$target_keep)
  keep_file_subset<-keep_file[sample(1000),]
  write.table(keep_file_subset, paste0(opt$output_dir,'_subset.keep'), col.names=F, row.names=F, quote=F)

}

rm(keep_file,keep_file_subset, fam_file)
gc()

if(is.na(opt$target_fam)){
  if(is.na(opt$target_plink)){
    for(i in 1:22){
      system(paste0(opt$plink,' --bfile ',opt$target_plink_chr,i,' --keep ',opt$output_dir,'_subset.keep --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --exclude ',opt$output_dir,'long_ld.exclude --indep-pairwise 1000 5 0.2 --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
    }
  } else {
    system(paste0(opt$plink,' --bfile ',opt$target_plink,' --keep ',opt$output_dir,'_subset.keep --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --exclude ',opt$output_dir,'long_ld.exclude --indep-pairwise 1000 5 0.2 --chr 1-22 --write-snplist --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
  }
} else {
  if(is.na(opt$target_plink)){
    for(i in 1:22){
      system(paste0(opt$plink,' --bfile ',opt$target_plink_chr,i,' --fam ',opt$target_fam,' --keep ',opt$output_dir,'_subset.keep --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --exclude ',opt$output_dir,'long_ld.exclude --indep-pairwise 1000 5 0.2 --write-snplist --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
    }
  } else {
    system(paste0(opt$plink,' --bfile ',opt$target_plink,' --fam ',opt$target_fam,' --keep ',opt$output_dir,'_subset.keep --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --exclude ',opt$output_dir,'long_ld.exclude --indep-pairwise 1000 5 0.2 --chr 1-22 --write-snplist --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
  }
}

# If there are over 100,000 markers remaining, restrict to 100,000
if(is.na(opt$target_plink)){
  prune_all<-NULL
  for(i in 1:22){
    prune_all<-rbind(prune_all, fread(paste0(opt$output_dir,'target.QC.chr',i,'.prune.in'), header=F))
    prune_all<-prune_all[sample(1:100000),]
  }
} else {
  prune_all<-fread(paste0(opt$output_dir,'target.QC.prune.in'), header=F)
  prune_all<-prune_all[sample(1:100000),]
}

write.table(prune_all, paste0(opt$output_dir,'target.QC.prune.in'), col.names=F, row.names=F, quote=F)

# Apply relatedness threshold
if(opt$save_grm == T){
  if(is.na(opt$target_keep)){
    if(is.na(opt$target_fam)){
      if(is.na(opt$target_plink)){
        # Create merge list
        ref_merge_list<-paste0(opt$target_plink,1:22)
        write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)
        
        system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --threads 1 --rel-cutoff ',opt$rel_thresh,' --make-grm-bin --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
      } else {
        system(paste0(opt$plink,' --bfile ',opt$target_plink,' --threads 1  --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --make-grm-bin --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
      }
    } else {
      if(is.na(opt$target_plink)){
        # Create merge list
        ref_merge_list<-paste0(opt$target_plink,1:22)
        write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)
        
        system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --fam ',opt$target_fam,' --threads 1  --rel-cutoff ',opt$rel_thresh,' --make-grm-bin --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
      } else {
        system(paste0(opt$plink,' --bfile ',opt$target_plink,' --fam ',opt$target_fam,' --threads 1  --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --make-grm-bin --chr 1-22 --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
      }
    }
  } else {
    if(is.na(opt$target_fam)){
      if(is.na(opt$target_plink)){
        # Create merge list
        ref_merge_list<-paste0(opt$target_plink,1:22)
        write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)
        
        system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --make-grm-bin --keep ',opt$target_keep,' --threads 1 --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
      } else {
        system(paste0(opt$plink,' --bfile ',opt$target_plink,' --make-grm-bin --keep ',opt$target_keep,' --threads 1  --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
      }
    } else {
      if(is.na(opt$target_plink)){
        # Create merge list
        ref_merge_list<-paste0(opt$target_plink,1:22)
        write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)
        
        system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --fam ',opt$target_fam,' --keep ',opt$target_keep,' --threads 1 --make-grm-bin --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
      } else {
        system(paste0(opt$plink,' --bfile ',opt$target_plink,' --fam ',opt$target_fam,' --keep ',opt$target_keep,' --threads 1  --rel-cutoff ',opt$rel_thresh,' --make-grm-bin --extract ',opt$output_dir,'target.QC.prune.in --chr 1-22 --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
      }
    }
  }
} else {
  if(is.na(opt$target_keep)){
    if(is.na(opt$target_fam)){
      if(is.na(opt$target_plink)){
        # Create merge list
        ref_merge_list<-paste0(opt$target_plink,1:22)
        write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)
        
        system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --threads 1 --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
      } else {
        system(paste0(opt$plink,' --bfile ',opt$target_plink,' --threads 1  --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
      }
    } else {
      if(is.na(opt$target_plink)){
        # Create merge list
        ref_merge_list<-paste0(opt$target_plink,1:22)
        write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)
        
        system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --fam ',opt$target_fam,' --threads 1  --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
      } else {
        system(paste0(opt$plink,' --bfile ',opt$target_plink,' --fam ',opt$target_fam,' --threads 1  --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --chr 1-22 --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
      }
    }
  } else {
    if(is.na(opt$target_fam)){
      if(is.na(opt$target_plink)){
        # Create merge list
        ref_merge_list<-paste0(opt$target_plink,1:22)
        write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)
        
        system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --keep ',opt$target_keep,' --threads 1 --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
      } else {
        system(paste0(opt$plink,' --bfile ',opt$target_plink,' --keep ',opt$target_keep,' --threads 1  --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
      }
    } else {
      if(is.na(opt$target_plink)){
        # Create merge list
        ref_merge_list<-paste0(opt$target_plink,1:22)
        write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)
        
        system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --fam ',opt$target_fam,' --keep ',opt$target_keep,' --threads 1  --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
      } else {
        system(paste0(opt$plink,' --bfile ',opt$target_plink,' --fam ',opt$target_fam,' --keep ',opt$target_keep,' --threads 1  --rel-cutoff ',opt$rel_thresh,' --extract ',opt$output_dir,'target.QC.prune.in --chr 1-22 --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
      }
    }
  }
}
###
# Report the number of individuals removed by relatedness threshold
###

unrelated<-fread(paste0(opt$output_dir,'target.QC.rel.id'))

write.table(unrelated, paste0(opt$output,'.unrelated.keep'), col.names=F, row.names=F, quote=F)

system(paste0('rm ',opt$output_dir, 'target.*'))
system(paste0('rm ',opt$output_dir, 'long_ld.exclude'))
system(paste0('rm ',opt$output_dir,'target.QC.rel.id'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
