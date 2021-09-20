#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--target_plink_chr", action="store", default=NA, type='character',
    help="Path to per chromosome target PLINK files [required]"),
make_option("--target_plink", action="store", default=NA, type='character',
    help="Path to per chromosome target PLINK files [required]"),
make_option("--ref_plink_chr", action="store", default=NA, type='character',
    help="Path to per chromosome reference PLINK files [required]"),
make_option("--target_fam", action="store", default=NA, type='character',
    help="Target sample fam file. [optional]"),
make_option("--maf", action="store", default=0.05, type='numeric',
    help="Minor allele frequency threshold [optional]"),
make_option("--geno", action="store", default=0.02, type='numeric',
    help="Variant missingness threshold [optional]"),
make_option("--hwe", action="store", default=1e-6, type='numeric',
    help="Hardy Weinberg p-value threshold. [optional]"),
make_option("--n_pcs", action="store", default=10, type='numeric',
		help="Number of PCs (min=4) [optional]"),
make_option("--plink", action="store", default='plink', type='character',
		help="Path PLINK software binary [required]"),
make_option("--plink2", action="store", default='plink', type='character',
		help="Path PLINK software binary [required]"),
make_option("--output", action="store", default='./PC_projector_output/Output', type='character',
		help="Path for output files [required]"),
make_option("--ref_pop_scale", action="store", default=NA, type='character',
		help="List of keep files ancestry specific scaling [optional]"),    
make_option("--pop_data", action="store", default=NA, type='character',
    help="Population data for the reference samples [optional]"),    
make_option("--model_method", action="store", default='glmnet', type='character',
    help="Method used for generate prediction model [optional]"),    
make_option("--SD_rule", action="store", default=F, type='logical',
    help="Logical indicating whether the 3SD rule should be used to define ancestry, or the model-based approach [optional]"),    
make_option("--prob_thresh", action="store", default='NA', type='numeric',
    help="Indicates whether probability threshold should be used when defining ancestry [optional]"),    
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(caret)
library(pROC)
library(verification)
library(ggplot2)
library(cowplot)

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# Ancestry_identifier.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

if(is.na(opt$ref_pop_scale)){
	sink(file = paste(opt$output,'.log',sep=''), append = F)
	cat('ref_pop must be specified\n')
	sink()
	q()
}

###########
# Perform QC of target sample genotypes
###########

if(is.na(opt$target_fam)){
  if(is.na(opt$target_plink)){
    for(i in 1:22){
      system(paste0(opt$plink,' --bfile ',opt$target_plink_chr,i,' --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --write-snplist --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
    }
  } else {
    system(paste0(opt$plink,' --bfile ',opt$target_plink,' --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --write-snplist --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
  }
} else {
  if(is.na(opt$target_plink)){
    for(i in 1:22){
      system(paste0(opt$plink,' --bfile ',opt$target_plink_chr,i,' --fam ',opt$target_fam,' --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --write-snplist --out ',opt$output_dir,'target.QC.chr',i,' --memory ',floor(opt$memory*0.7)))
    }
  } else {
      system(paste0(opt$plink,' --bfile ',opt$target_plink,' --fam ',opt$target_fam,' --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --write-snplist --out ',opt$output_dir,'target.QC --memory ',floor(opt$memory*0.7)))
  }
}

if(is.na(opt$target_plink)){
  target_qc_snplist<-NULL
  for(i in 1:22){
    target_qc_snplist_i<-fread(paste0(opt$output_dir,'target.QC.chr',i,'.snplist'), header=F)
    target_qc_snplist<-rbind(target_qc_snplist, target_qc_snplist_i)
  }
  
  write.table(target_qc_snplist$V1, paste0(opt$output_dir,'target_QC.snplist'), row.names=F, col.names=F, quote=F)
  
  system(paste0('rm ',paste0(opt$output_dir,'target.QC.chr*')))
  
} else {
  target_qc_snplist<-fread(paste0(opt$output_dir,'target.QC.snplist'), header=F)

  write.table(target_qc_snplist$V1, paste0(opt$output_dir,'target_QC.snplist'), row.names=F, col.names=F, quote=F)
  
  system(paste0('rm ',paste0(opt$output_dir,'target.QC.*')))
}

###
# Merge the per chromosome reference genetic data, retaining only comp_snp
###

# Extract intersect with target from the reference files which also meet QC
for(i in 1:22){
  system(paste0(opt$plink,' --bfile ', opt$ref_plink_chr,i,' --threads 1 --geno ',opt$geno,' --maf ',opt$maf,' --hwe ',opt$hwe,' --make-bed --extract ',opt$output_dir,'target_QC.snplist --out ',opt$output_dir,'ref_intersect_chr',i,' --memory ',floor(opt$memory*0.7)))
}

##
# Check alleles match between reference and target
##

# This is important in instances when the target and reference have not already been harmonised
ref_bim<-NULL
for(i in 1:22){
  ref_bim<-rbind(ref_bim, fread(paste0(opt$output_dir,'ref_intersect_chr',i,'.bim')))
}

ref_bim<-ref_bim[,c('V1','V2','V4','V5','V6')]
names(ref_bim)<-c('CHR','SNP','BP','A1','A2')

if(is.na(opt$target_plink)){
  targ_bim<-NULL
  for(i in 1:22){
    targ_bim<-rbind(targ_bim, fread(paste0(opt$target_plink_chr,i,'.bim')))
  }
} else {
  targ_bim<-fread(paste0(opt$target_plink,'.bim'))
}

targ_bim<-targ_bim[,c('V1','V2','V4','V5','V6')]
names(targ_bim)<-c('CHR','SNP','BP','A1','A2')

# Create IUPAC codes in target data
targ_bim$IUPAC[targ_bim$A1 == 'A' & targ_bim$A2 =='T' | targ_bim$A1 == 'T' & targ_bim$A2 =='A']<-'W'
targ_bim$IUPAC[targ_bim$A1 == 'C' & targ_bim$A2 =='G' | targ_bim$A1 == 'G' & targ_bim$A2 =='C']<-'S'
targ_bim$IUPAC[targ_bim$A1 == 'A' & targ_bim$A2 =='G' | targ_bim$A1 == 'G' & targ_bim$A2 =='A']<-'R'
targ_bim$IUPAC[targ_bim$A1 == 'C' & targ_bim$A2 =='T' | targ_bim$A1 == 'T' & targ_bim$A2 =='C']<-'Y'
targ_bim$IUPAC[targ_bim$A1 == 'G' & targ_bim$A2 =='T' | targ_bim$A1 == 'T' & targ_bim$A2 =='G']<-'K'
targ_bim$IUPAC[targ_bim$A1 == 'A' & targ_bim$A2 =='C' | targ_bim$A1 == 'C' & targ_bim$A2 =='A']<-'M'
targ_bim$SNP_IUPAC<-paste0(targ_bim$SNP,':',targ_bim$IUPAC)

# Create IUPAC codes in ref data
ref_bim$IUPAC[ref_bim$A1 == 'A' & ref_bim$A2 =='T' | ref_bim$A1 == 'T' & ref_bim$A2 =='A']<-'W'
ref_bim$IUPAC[ref_bim$A1 == 'C' & ref_bim$A2 =='G' | ref_bim$A1 == 'G' & ref_bim$A2 =='C']<-'S'
ref_bim$IUPAC[ref_bim$A1 == 'A' & ref_bim$A2 =='G' | ref_bim$A1 == 'G' & ref_bim$A2 =='A']<-'R'
ref_bim$IUPAC[ref_bim$A1 == 'C' & ref_bim$A2 =='T' | ref_bim$A1 == 'T' & ref_bim$A2 =='C']<-'Y'
ref_bim$IUPAC[ref_bim$A1 == 'G' & ref_bim$A2 =='T' | ref_bim$A1 == 'T' & ref_bim$A2 =='G']<-'K'
ref_bim$IUPAC[ref_bim$A1 == 'A' & ref_bim$A2 =='C' | ref_bim$A1 == 'C' & ref_bim$A2 =='A']<-'M'
ref_bim$SNP_IUPAC<-paste0(ref_bim$SNP,':',ref_bim$IUPAC)

# Merge target and reference based on SNP id
ref_target<-merge(ref_bim, targ_bim, by='SNP')

# Identify SNPs for which alleles need to be flipped
flip_tmp<-ref_target[(ref_target$IUPAC.x == 'R' & ref_target$IUPAC.y == 'Y' | 
                        ref_target$IUPAC.x == 'Y' & ref_target$IUPAC.y == 'R' | 
                        ref_target$IUPAC.x == 'K' & ref_target$IUPAC.y == 'M' |
                        ref_target$IUPAC.x == 'M' & ref_target$IUPAC.y == 'K'),]

# Idenitfy SNPs which match the reference alleles
incl<-ref_target[ ref_target$IUPAC.x == 'R' & ref_target$IUPAC.y == 'R' | 
                    ref_target$IUPAC.x == 'Y' & ref_target$IUPAC.y == 'Y' | 
                    ref_target$IUPAC.x == 'K' & ref_target$IUPAC.y == 'K' |
                    ref_target$IUPAC.x == 'M' & ref_target$IUPAC.y == 'M' ]

# If a SNP that needs to be flipped has a duplicate that is on the correct strand, remove it.
flip<-flip_tmp[!(flip_tmp$SNP %in% incl$SNP)]

# Combine SNPs that match and those that need to be flipped.
incl<-rbind(incl,flip)

if(dim(flip)[1] > 0){
  write.table(flip$SNP, paste0(opt$output_dir,'ref_flip.snplist'), col.names=F, row.names=F, quote=F)
}

write.table(incl$SNP, paste0(opt$output_dir,'ref_allele_match.snplist'), col.names=F, row.names=F, quote=F)

# Merge subset reference
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Merging per chromosome reference data...')
sink()

# Create merge list
ref_merge_list<-paste0(opt$output_dir,'ref_intersect_chr',1:22)

write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)

# Merge
if(dim(flip)[1] > 0){
  system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --extract ',opt$output_dir,'ref_allele_match.snplist --flip ',opt$output_dir,'ref_flip.snplist --threads 1 --make-bed --out ',opt$output_dir,'ref_merge --memory ',floor(opt$memory*0.7)))
} else {
  system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --extract ',opt$output_dir,'ref_allele_match.snplist --threads 1 --make-bed --out ',opt$output_dir,'ref_merge --memory ',floor(opt$memory*0.7)))
}
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

# Delete temporary per chromosome reference files
system(paste0('rm ',opt$output_dir,'ref_intersect_chr*'))
if(dim(flip)[1] > 0){
  system(paste0('rm ',opt$output_dir,'ref_flip.snplist'))
}
system(paste0('rm ',opt$output_dir,'ref_allele_match.snplist'))

###
# Create SNP list for LD pruning
# Remove regions of long range LD which can confound estimates of ancestry estimates (REF: PMC2443852)
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Identifying LD independent SNPs based on reference data...')
sink()

# Read in the bim file
ref_bim<-data.frame(fread(paste0(opt$output_dir,'ref_merge.bim')))

# Create file removing these regions.
long_ld_exclude<-ref_bim$V2[ (ref_bim$V1 == 1 & ref_bim$V4 >= 48e6 & ref_bim$V4 <= 52e6) |
                                  (ref_bim$V1 == 2 & ref_bim$V4 >= 86e6 & ref_bim$V4 <= 100.5e6) |
                                  (ref_bim$V1 == 2 & ref_bim$V4 >= 134.5e6 & ref_bim$V4 <= 138e6) |
                                  (ref_bim$V1 == 2 & ref_bim$V4 >= 183e6 & ref_bim$V4 <= 190e6) |
                                  (ref_bim$V1 == 3 & ref_bim$V4 >= 47.5e6 & ref_bim$V4 <= 50e6) |
                                  (ref_bim$V1 == 3 & ref_bim$V4 >= 83.5e6 & ref_bim$V4 <= 87e6) |
                                  (ref_bim$V1 == 3 & ref_bim$V4 >= 89e6 & ref_bim$V4 <= 97.5e6) |
                                  (ref_bim$V1 == 5 & ref_bim$V4 >= 44.5e6 & ref_bim$V4 <= 50.5e6) |
                                  (ref_bim$V1 == 5 & ref_bim$V4 >= 98e6 & ref_bim$V4 <= 100.5e6) |
                                  (ref_bim$V1 == 5 & ref_bim$V4 >= 129e6 & ref_bim$V4 <= 132e6) |
                                  (ref_bim$V1 == 5 & ref_bim$V4 >= 135.5e6 & ref_bim$V4 <= 138.5e6) |
                                  (ref_bim$V1 == 6 & ref_bim$V4 >= 25.5e6 & ref_bim$V4 <= 33.5e6) |
                                  (ref_bim$V1 == 6 & ref_bim$V4 >= 57e6 & ref_bim$V4 <= 64e6) |
                                  (ref_bim$V1 == 6 & ref_bim$V4 >= 140e6 & ref_bim$V4 <= 142.5e6) |
                                  (ref_bim$V1 == 7 & ref_bim$V4 >= 55e6 & ref_bim$V4 <= 66e6) |
                                  (ref_bim$V1 == 8 & ref_bim$V4 >= 8e6 & ref_bim$V4 <= 12e6) |
                                  (ref_bim$V1 == 8 & ref_bim$V4 >= 43e6 & ref_bim$V4 <= 50e6) |
                                  (ref_bim$V1 == 8 & ref_bim$V4 >= 112e6 & ref_bim$V4 <= 115e6) |
                                  (ref_bim$V1 == 10 & ref_bim$V4 >= 37e6 & ref_bim$V4 <= 43e6) |
                                  (ref_bim$V1 == 11 & ref_bim$V4 >= 46e6 & ref_bim$V4 <= 57e6) |
                                  (ref_bim$V1 == 11 & ref_bim$V4 >= 87.5e6 & ref_bim$V4 <= 90.5e6) |
                                  (ref_bim$V1 == 12 & ref_bim$V4 >= 33e6 & ref_bim$V4 <= 40e6) |
                                  (ref_bim$V1 == 12 & ref_bim$V4 >= 109.5e6 & ref_bim$V4 <= 112e6) |
                                  (ref_bim$V1 == 20 & ref_bim$V4 >= 32e6 & ref_bim$V4 <= 34.5e6)]

write.table(long_ld_exclude, paste0(opt$output_dir,'long_ld.exclude'), col.names=F, row.names=F, quote=F)
  
# Identify LD independent SNPs.
system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_merge --threads 1 --exclude ',opt$output_dir,'long_ld.exclude --indep-pairwise 1000 5 0.2 --out ',opt$output_dir,'ref_merge --memory ',floor(opt$memory*0.7)))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Calculate PCs in the reference sample for scaling the target sample factor scores.
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Computing reference PCs...')
sink()

# Extract LD independent SNPs
system(paste0(opt$plink,' --bfile ',opt$output_dir,'ref_merge --threads 1 --extract ',opt$output_dir,'ref_merge.prune.in --make-bed --out ',opt$output_dir,'ref_merge_pruned --memory ',floor(opt$memory*0.7)))

# Calculate SNP weights
system(paste0(opt$plink2,' --bfile ',opt$output_dir,'ref_merge_pruned --threads 1 --pca ',opt$n_pcs,' biallelic-var-wts  --out ',opt$output_dir,'ref_merge_pruned --memory ',floor(opt$memory*0.7)))

# Calculate PCs in the reference
system(paste0(opt$plink2,' --bfile ',opt$output_dir,'ref_merge_pruned --threads 1 --score ',opt$output_dir,'ref_merge_pruned.eigenvec.var header-read 2 3 --score-col-nums 5-',as.numeric(opt$n_pcs)+4,' --out ',opt$output_dir,'ref_merge_pruned_score --memory ',floor(opt$memory*0.7)))

# Read in reference PC scores
PCs_ref<-data.frame(fread(paste0(opt$output_dir,'ref_merge_pruned_score.sscore')))
PCs_ref<-PCs_ref[,c(1:2,5:dim(PCs_ref)[2])]
names(PCs_ref)<-c('FID','IID',paste0('PC',1:as.numeric(opt$n_pcs)))

fwrite(PCs_ref, paste0(opt$output,'.eigenvec'), sep=' ')

# Scale across all individuals
PCs_ref_centre_scale<-data.frame(PC=names(PCs_ref[-1:-2]),
							  Mean=sapply(PCs_ref[,-1:-2], function(x) mean(x)),
							  SD=sapply(PCs_ref[,-1:-2], function(x) sd(x)),
							  row.names=seq(1:as.numeric(opt$n_pcs)))

fwrite(PCs_ref_centre_scale, paste0(opt$output,'.scale'), sep=' ')

rm(PCs_ref_centre_scale)
gc()

if(!is.na(opt$ref_pop_scale)){
  # Calculate the mean and sd of scores for each population specified in pop_scale
  pop_keep_files<-read.table(opt$ref_pop_scale, header=F, stringsAsFactors=F)

  for(k in 1:dim(pop_keep_files)[1]){
  	pop<-pop_keep_files$V1[k]
  	keep<-fread(pop_keep_files$V2[k], header=F)
  	PCs_ref_keep<-PCs_ref[(PCs_ref$FID %in% keep$V1),]

    PCs_ref_centre_scale<-data.frame(PC=names(PCs_ref_keep[-1:-2]),
    								  Mean=sapply(PCs_ref_keep[,-1:-2], function(x) mean(x)),
    								  SD=sapply(PCs_ref_keep[,-1:-2], function(x) sd(x)),
    								  row.names=seq(1:opt$n_pcs))

  	fwrite(PCs_ref_centre_scale, paste0(opt$output,'.',pop,'.scale'), sep=' ')
	
	rm(PCs_ref_centre_scale)
	gc()
	}
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Create model predicting ref_pop groups
###

if(!is.na(opt$ref_pop_scale)){

	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Deriving model predicting ref_pop groups...')
	sink()

	# Read in whole sample scale file
	PCs_ref_centre_scale<-fread(paste0(opt$output,'.scale'))

	# Scale the reference PCs
	PCs_ref_scaled<-PCs_ref
	for(i in 1:dim(PCs_ref_centre_scale)[1]){
		PCs_ref_scaled[[paste0('PC',i)]]<-PCs_ref[[paste0('PC',i)]]-PCs_ref_centre_scale$Mean[PCs_ref_centre_scale$PC == paste0('PC',i)]
		PCs_ref_scaled[[paste0('PC',i)]]<-PCs_ref_scaled[[paste0('PC',i)]]/PCs_ref_centre_scale$SD[PCs_ref_centre_scale$PC == paste0('PC',i)]
	}
	
	# Label individuals with ref_pop groups
	pop<-NULL
	for(i in 1:dim(pop_keep_files)[1]){
		keep<-fread(pop_keep_files$V2[i], header=F)
		keep$pop<-pop_keep_files$V1[i]
		pop<-rbind(pop,keep)
	}
	names(pop)<-c('FID','IID','pop')
	PCs_ref_scaled_pop<-merge(PCs_ref_scaled,pop, by=c('FID','IID'))
	rm(PCs_ref_scaled)
	gc()
	
	# Build model
	model <- train(y=as.factor(PCs_ref_scaled_pop$pop), x=PCs_ref_scaled_pop[grepl('PC',names(PCs_ref_scaled_pop))], method=opt$model_method, metric='logLoss', trControl=trainControl(method="cv", number=5, classProbs= TRUE, savePredictions = 'final', summaryFunction = multiClassSummary))
	
	# Save performance information
	sink(file = paste(opt$output,'.pop_model_prediction_details.txt',sep=''), append = F)
	print(model)
	cat('\n')
	obs_pre_tab<-table(model$pred$obs, model$pred$pred)
	dimnames(obs_pre_tab)<-list(paste('obs',dimnames(obs_pre_tab)[[1]]),paste('pred',dimnames(obs_pre_tab)[[2]]))
	
	# Show confusion matrix before and after applying probability threshold
	cat('Confusion matrix without threshold:\n')
	print(obs_pre_tab)
	
	if(!is.na(opt$prob_thresh)){
	  model$pred$max_prob<-apply(model$pred[,unique(PCs_ref_scaled_pop$pop)], 1, max)
	  model$pred<-model$pred[model$pred$max_prob > opt$prob_thresh,]
	  
	  obs_pre_tab_thresh<-table(model$pred$obs, model$pred$pred)
	  dimnames(obs_pre_tab_thresh)<-list(paste('obs',dimnames(obs_pre_tab_thresh)[[1]]),paste('pred',dimnames(obs_pre_tab_thresh)[[2]]))
	  
  	cat('\n')
  	cat(paste0('Confusion matrix with ',opt$prob_thresh,' threshold:\n'))
  	print(obs_pre_tab_thresh)
	}
	
	sink()
	
	saveRDS(model$finalModel, paste0(opt$output,'.pop_model.rds'))
	
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Done!\n')
	sink()
}

###
# Rename files
###

system(paste0('mv ',opt$output_dir,'ref_merge_pruned.eigenvec.var ',opt$output,'.eigenvec.var'))
system(paste0('rm ',opt$output_dir,'ref_merge*'))
system(paste0('rm ',opt$output_dir,'long_ld.exclude'))
system(paste0('rm ',opt$output_dir,'target_QC.snplist'))

#####
# Calculate PCs in target sample
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating PCs in the target sample...')
sink()

system(paste0('cut -f 2 ',opt$output,'.eigenvec.var | tail -n +2 > ',opt$output_dir,'score_file.snplist'))

if(is.na(opt$target_fam)){
  if(is.na(opt$target_plink)){
    for(i in 1:22){
      system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --extract ',opt$output_dir,'score_file.snplist --score ',opt$output,'.eigenvec.var header-read 2 3 no-mean-imputation --threads 1 --score-col-nums 5-',as.numeric(opt$n_pcs)+4,' --out ',opt$output_dir,'profiles.chr',i,' --memory ',floor(opt$memory*0.9)))
    }
  } else {
    system(paste0(opt$plink2, ' --bfile ',opt$target_plink,' --extract ',opt$output_dir,'score_file.snplist --score ',opt$output,'.eigenvec.var header-read 2 3 no-mean-imputation --threads 1 --score-col-nums 5-',as.numeric(opt$n_pcs)+4,' --out ',opt$output_dir,'profiles --memory ',floor(opt$memory*0.9)))
  }
} else {
  if(is.na(opt$target_plink)){
    for(i in 1:22){
      system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --fam ',opt$target_fam,' --extract ',opt$output_dir,'score_file.snplist --score ',opt$output,'.eigenvec.var header-read 2 3 no-mean-imputation --threads 1 --score-col-nums 5-',as.numeric(opt$n_pcs)+4,' --out ',opt$output_dir,'profiles.chr',i,' --memory ',floor(opt$memory*0.9)))
    }
  } else {
    system(paste0(opt$plink2, ' --bfile ',opt$target_plink,' --fam ',opt$target_fam,' --extract ',opt$output_dir,'score_file.snplist --score ',opt$output,'.eigenvec.var header-read 2 3 no-mean-imputation --threads 1 --score-col-nums 5-',as.numeric(opt$n_pcs)+4,' --out ',opt$output_dir,'profiles --memory ',floor(opt$memory*0.9)))
  }
}

system(paste0('rm ',opt$output_dir,'score_file.snplist'))

if(is.na(opt$target_plink)){
  # Add up the scores across chromosomes
  scores<-fread(cmd=paste0('cut -f 1-2 ',opt$output_dir,'profiles.chr22.sscore'))
  names(scores)<-c('FID','IID')
  
  var_list<-fread(paste0(opt$output,'.eigenvec.var'))
  nsnp_all<-0
  for(i in 1:22){
  	profile<-data.frame(fread(paste0(opt$output_dir,'profiles.chr',i,'.sscore')))
  	profile<-as.matrix(profile[,grepl('PC',names(profile))])
  	bim<-fread(paste0(opt$target_plink_chr,i,'.bim'))
  	nsnp<-sum(bim$V2 %in% var_list$ID)
  	nsnp_all<-nsnp_all+nsnp
  	profile<-profile*nsnp
  	if(i == 1){
  		profile_all<-profile
  	} else {
  		profile_all<-profile_all+profile	
  	}
  print(i)
  }
  
  profile_all<-profile_all/nsnp_all
} else {
  # Add up the scores across chromosomes
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
rm(scores,profile_all,var_list,nsnp_all)
gc()

###
# Clean up temporary files
###

system(paste0('rm ',opt$output_dir,'profiles*'))

###
# Create plot PC scores of target sample compared to the reference
###
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Plotting target sample PCs on reference...')
sink()

# Read in population data
pop_data<-data.frame(fread(opt$pop_data))
names(pop_data)[1]<-'IID'
pop_data$FID<-pop_data$IID

# Read in reference sample PCs
ref_PCs<-data.frame(fread(paste0(opt$output,'.eigenvec')))
ref_PCs<-merge(ref_PCs, pop_data, by=c('FID','IID'))

# Insert pop_data columns into target PCs
new_cols<-names(ref_PCs[!grepl('PC|ID', names(ref_PCs))])
new_cols_2<-data.frame(matrix(rep('Target',length(new_cols)),ncol=length(new_cols)))
names(new_cols_2)<-names(ref_PCs[!grepl('PC|ID', names(ref_PCs))])
targ_PCs<-cbind(targ_PCs,new_cols_2)

# Combine the two sets
ref_PCs_targ_PCs<-rbind(ref_PCs,targ_PCs)

rm(ref_PCs)
gc()

Label_groups<-names(ref_PCs_targ_PCs[!grepl('PC|IID|FID',names(ref_PCs_targ_PCs))])

for(i in Label_groups){
PC_1_2<-ggplot(ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] != 'Target',], aes(x=PC1,y=PC2, colour=get(i))) + 
  geom_point() + 
	geom_point(data=ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] == 'Target',], aes(x=PC1,y=PC2), colour='black', shape=21) + 
  ggtitle("PCs 1 and 2") +
	labs(colour="")
PC_3_4<-ggplot(ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] != 'Target',], aes(x=PC3,y=PC4, colour=get(i))) + 
  geom_point() + 
	geom_point(data=ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] == 'Target',], aes(x=PC3,y=PC4), colour='black', shape=21) + 
  ggtitle("PCs 3 and 4") +
	labs(colour="")


png(paste0(opt$output,'.PCs_plot_',i,'.png'), units='px', res=300, width=4000, height=2500)
print(plot_grid(PC_1_2,PC_3_4))
dev.off()

rm(PC_1_2,PC_3_4)
gc()

print(i)
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Estimate probability of outcomes in model
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Estimating probability of each population...')
sink()

# Read in the reference scale file
pop_model_scale<-fread(paste0(opt$output,'.scale'))

# Scale the target PCs
targ_PCs_scaled<-targ_PCs
for(i in 1:dim(pop_model_scale)[1]){
	targ_PCs_scaled[[paste0('PC',i)]]<-targ_PCs[[paste0('PC',i)]]-pop_model_scale$Mean[pop_model_scale$PC == paste0('PC',i)]
	targ_PCs_scaled[[paste0('PC',i)]]<-targ_PCs_scaled[[paste0('PC',i)]]/pop_model_scale$SD[pop_model_scale$PC == paste0('PC',i)]
}

# Read in model
pop_model<-readRDS(paste0(opt$output,'.pop_model.rds'))
pop_model_pred<-predict(object = pop_model, newx = data.matrix(targ_PCs_scaled[grepl('PC',names(targ_PCs_scaled))]), type = "response", s=pop_model$lambdaOpt)
pop_model_pred<-as.data.frame.table(pop_model_pred)
pop_model_pred<-data.table(	FID=targ_PCs_scaled$FID,
														IID=targ_PCs_scaled$IID,
														pop=as.character(pop_model_pred$Var2),
														prob=round(pop_model_pred$Freq,3))
		
pop_model_pred<-dcast.data.table(pop_model_pred, formula=FID + IID~pop, value.var = "prob")

fwrite(pop_model_pred, paste0(opt$output,'.model_pred'), sep='\t')

# Create keep files based on the results
if(!is.na(opt$prob_thresh)){
  pop_model_pred$max_prob<-apply(pop_model_pred[,-1:-2], 1, max)
  pop_model_pred<-pop_model_pred[pop_model_pred$max_prob > opt$prob_thresh,]
  pop_model_pred$max_prob<-NULL
}

N_group<-NULL
for(i in names(pop_model_pred[,-1:-2])){
	tmp_keep<-pop_model_pred[apply(pop_model_pred[,-1:-2], 1, function(x) x[i] == max(x)),1:2]
	N_group<-rbind(N_group, data.frame(Group=i, N=dim(tmp_keep)[1]))
	fwrite(tmp_keep, paste0(opt$output,'.model_pred.',i,'.keep'), sep=' ', col.names=F)
}

rm(targ_PCs_scaled,pop_model_pred)
gc()

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('----------\n')
cat('N per group based on model:\n')
print(N_group)
cat('----------\n')
sink()

###
# Identify individuals that are within 3SD of population specific mean for all PCs and write out scaled PCs.
###

targ_PCs<-targ_PCs[,grepl('FID|IID|PC',names(targ_PCs))]
# Read in pop_scale_for_keep
pop_scale_for_keep<-paste0(opt$output,'.',pop_keep_files$V1,'.scale')

N_group<-NULL
for(i in 1:length(pop_scale_for_keep)){
	# Idenitfy name of population based on scale file
	pop_name<-gsub('.scale','',substr(pop_scale_for_keep[i], nchar(pop_scale_for_keep[i])-9+1, nchar(pop_scale_for_keep[i])))
	pop_scale_for_keep_i<-fread(pop_scale_for_keep[i])
	
	# Scale the target based on population scale file
	targ_PCs_scaled_i<-targ_PCs[,grepl('FID|IID|PC', names(targ_PCs))]
	for(j in 1:dim(pop_scale_for_keep_i)[1]){
		targ_PCs_scaled_i[[paste0('PC',j)]]<-targ_PCs[[paste0('PC',j)]]-pop_scale_for_keep_i$Mean[pop_scale_for_keep_i$PC == paste0('PC',j)]
		targ_PCs_scaled_i[[paste0('PC',j)]]<-targ_PCs_scaled_i[[paste0('PC',j)]]/pop_scale_for_keep_i$SD[pop_scale_for_keep_i$PC == paste0('PC',j)]
		targ_PCs_scaled_i[[paste0('PC',j)]]<-round(targ_PCs_scaled_i[[paste0('PC',j)]],3)
	}
	
	# Remove anyone with a PC value >3 or -3 (i.e. 3SD from the population mean
	targ_PCs_scaled_i<-targ_PCs_scaled_i[!apply(targ_PCs_scaled_i[,-1:-2], 1, function(x) any(x > 3 | x < -3)),]
	
	N_group<-rbind(N_group, data.frame(Group=pop_name, N=dim(targ_PCs_scaled_i)[1]))
	
	# Save keep file of individuals that fit the population
	fwrite(targ_PCs_scaled_i[,1:2], paste0(opt$output,'.',pop_name,'.keep'), col.names=F, sep='\t')
	
	# Write the scaled PCs
	fwrite(targ_PCs_scaled_i, paste0(opt$output,'.',pop_name,'.eigenvec'), sep='\t')
	
	rm(pop_name,pop_scale_for_keep_i,targ_PCs_scaled_i)
	gc()
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('----------\n')
cat('N per group based on 3SD rule:\n')
print(N_group)
cat('----------\n')
sink()

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()

