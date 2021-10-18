#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--ref_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK files [required]"),
make_option("--ref_keep", action="store", default=NA, type='character',
		help="Keep file to subset individuals in reference for clumping [required]"),
make_option("--n_pcs", action="store", default=10, type='numeric',
		help="Number of PCs [optional]"),
make_option("--plink", action="store", default='plink', type='character',
		help="Path PLINK software binary [required]"),
make_option("--plink2", action="store", default='plink', type='character',
		help="Path PLINK2 software binary [required]"),
make_option("--output", action="store", default='./PC_projector_output/Output', type='character',
		help="Path for output files [required]"),
make_option("--ref_pop_scale", action="store", default=NA, type='character',
		help="List of keep files for scaling PCs [optional]"),    
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(caret)
library(pROC)
library(verification)

opt$output_dir<-dirname(opt$output)
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# ancestry_score_file_creator.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

###
# Merge the per chromosome reference genetic data and update IDs to be clearly distinguishable from the target
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Merging per chromosome reference data...')
sink()

# Create merge list
ref_merge_list<-paste0(opt$ref_plink_chr,1:22)

write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)

# Merge
if(!is.na(opt$ref_keep)){
  system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --threads 1 --make-bed --keep ',opt$ref_keep,' --out ',opt$output_dir,'ref_merge --memory ',floor(opt$memory*0.7)))
} else {
  system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --threads 1 --make-bed --out ',opt$output_dir,'ref_merge --memory ',floor(opt$memory*0.7)))  
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

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
system(paste0(opt$plink2,' --bfile ',opt$output_dir,'ref_merge_pruned --threads 1 --pca ',opt$n_pcs,' var-wts  --out ',opt$output_dir,'ref_merge_pruned --memory ',floor(opt$memory*0.7)))

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

if(!is.na(opt$ref_pop)){
  # Calculate the mean and sd of scores for each population specified in pop_scale
  pop_keep_files<-read.table(opt$ref_pop, header=F, stringsAsFactors=F)

  for(k in 1:dim(pop_keep_files)[1]){
  	pop<-pop_keep_files$V1[k]
  	keep<-fread(pop_keep_files$V2[k], header=F)
  	PCs_ref_keep<-PCs_ref[(PCs_ref$FID %in% keep$V1),]

    PCs_ref_centre_scale<-data.frame(PC=names(PCs_ref_keep[-1:-2]),
    								  Mean=sapply(PCs_ref_keep[,-1:-2], function(x) mean(x)),
    								  SD=sapply(PCs_ref_keep[,-1:-2], function(x) sd(x)),
    								  row.names=seq(1:100))

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

if(!is.na(opt$ref_pop)){

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
	
	# Build enet model
	enet_model <- train(y=as.factor(PCs_ref_scaled_pop$pop), x=PCs_ref_scaled_pop[grepl('PC',names(PCs_ref_scaled_pop))], method="glmnet", metric='logLoss', trControl=trainControl(method="cv", number=5, classProbs= TRUE, savePredictions = 'final', summaryFunction = multiClassSummary),tuneGrid = expand.grid(alpha = 0,lambda = 0))
	
	# Calculate the percentage of correctly individuals to each group		
	enet_model_n_correct<-NULL
	for(k in as.character(unique(enet_model$pred$obs))){
		tmp<-enet_model$pred[enet_model$pred$obs == k,]
	
	n_correct_tmp<-data.frame(	Group=k,
								N_obs=sum(tmp$obs == k),
								prop_correct=round(sum(tmp$obs == k & tmp$pred == k)/sum(tmp$obs == k),3))
	
	enet_model_n_correct<-rbind(enet_model_n_correct,n_correct_tmp)
	}
	
	write.table(enet_model_n_correct, paste0(opt$output,'.pop_enet_prediction_details.txt'), col.names=T, row.names=F, quote=F)

	saveRDS(enet_model$finalModel, paste0(opt$output,'.pop_enet_model.rds'))
	
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

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
