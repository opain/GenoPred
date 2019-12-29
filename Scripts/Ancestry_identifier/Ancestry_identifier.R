#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--target_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome target PLINK files [required]"),
make_option("--ref_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK files [required]"),
make_option("--n_pcs", action="store", default=10, type='numeric',
		help="Number of PCs [optional]"),
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
# Identify list of SNPs genotyped in all individuals in the target sample
###########

for(i in 1:22){
	system(paste0(opt$plink,' --bfile ',opt$target_plink_chr,i,' --threads 1 --missing --out ',opt$output_dir,'target.chr',i,' --memory ',floor(opt$memory*0.7)))
}

lmiss_all<-NULL
for(i in 1:22){
	lmiss<-fread(paste0(opt$output_dir,'target.chr',i,'.lmiss'))
	lmiss_all<-rbind(lmiss_all, lmiss)
}

comp_snp<-lmiss_all$SNP[lmiss_all$F_MISS < 0.02]

write.table(comp_snp, paste0(opt$output_dir,'target_comp.snplist'), row.names=F, col.names=F, quote=F)

###
# Merge the per chromosome reference genetic data, retaining only comp_snp
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Merging per chromosome reference data...')
sink()

# Create merge list
ref_merge_list<-paste0(opt$ref_plink_chr,1:22)

write.table(ref_merge_list, paste0(opt$output_dir,'ref_mergelist.txt'), row.names=F, col.names=F, quote=F)

# Merge
system(paste0(opt$plink,' --merge-list ',opt$output_dir,'ref_mergelist.txt --threads 1 --make-bed --extract ',opt$output_dir,'target_comp.snplist --out ',opt$output_dir,'ref_merge --memory ',floor(opt$memory*0.7)))
  
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
system(paste0('rm ',opt$output_dir,'target.chr*.lmiss'))
system(paste0('rm ',opt$output_dir,'target.chr*.imiss'))
system(paste0('rm ',opt$output_dir,'target.chr*.log'))
system(paste0('rm ',opt$output_dir,'target_comp.snplist'))

#####
# Calculate PCs in target sample
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating PCs in the target sample...')
sink()

system(paste0('cut -f 2 ',opt$output,'.eigenvec.var | tail -n +2 > ',opt$output_dir,'score_file.snplist'))

for(i in 1:22){
	system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --extract ',opt$output_dir,'score_file.snplist --score ',opt$output,'.eigenvec.var header-read 2 3 no-mean-imputation --threads 1 --score-col-nums 5-',as.numeric(opt$n_pcs)+4,' --out ',opt$output_dir,'profiles.chr',i,' --memory ',floor(opt$memory*0.9)))
}

system(paste0('rm ',opt$output_dir,'score_file.snplist'))

# Add up the scores across chromosomes
scores<-fread(paste0('cut -f 1-2 ',opt$output_dir,'profiles.chr22.sscore'))
names(scores)<-c('FID','IID')

var_list<-fread(paste0(opt$output,'.eigenvec.var'))
nsnp_all<-0
for(i in 1:22){
	profile<-data.frame(fread(paste0(opt$output_dir,'profiles.chr',i,'.sscore')))
	profile<-as.matrix(profile[,5:dim(profile)[2]])
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
PC_5_6<-ggplot(ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] != 'Target',], aes(x=PC5,y=PC6, colour=get(i))) + 
  geom_point() + 
	geom_point(data=ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] == 'Target',], aes(x=PC5,y=PC6), colour='black', shape=21) + 
  ggtitle("PCs 5 and 6") +
	labs(colour="")
PC_7_8<-ggplot(ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] != 'Target',], aes(x=PC7,y=PC8, colour=get(i))) + 
  geom_point() + 
	geom_point(data=ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] == 'Target',], aes(x=PC7,y=PC8), colour='black', shape=21) + 
  ggtitle("PCs 7 and 8") +
	labs(colour="")
PC_9_10<-ggplot(ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] != 'Target',], aes(x=PC9,y=PC10, colour=get(i))) + 
  geom_point() + 
	geom_point(data=ref_PCs_targ_PCs[ref_PCs_targ_PCs[[i]] == 'Target',], aes(x=PC9,y=PC10), colour='black', shape=21) + 
  ggtitle("PCs 9 and 10") +
	labs(colour="")

png(paste0(opt$output,'.PCs_plot_',i,'.png'), units='px', res=300, width=4000, height=2500)
print(plot_grid(PC_1_2,PC_3_4,PC_5_6,PC_7_8,PC_9_10))
dev.off()

rm(PC_1_2,PC_3_4,PC_5_6,PC_7_8,PC_9_10)
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
pop_model<-readRDS(paste0(opt$output,'.pop_enet_model.rds'))
pop_model_pred<-predict(object = pop_model, newx = data.matrix(targ_PCs_scaled[grepl('PC',names(targ_PCs_scaled))]), type = "response", s=pop_model$lambdaOpt)
pop_model_pred<-as.data.frame.table(pop_model_pred)
pop_model_pred<-data.table(	FID=targ_PCs_scaled$FID,
														IID=targ_PCs_scaled$IID,
														pop=as.character(pop_model_pred$Var2),
														prob=round(pop_model_pred$Freq,3))
							
pop_model_pred<-dcast.data.table(pop_model_pred, formula=FID + IID~pop, value.var = "prob")

fwrite(pop_model_pred, paste0(opt$output,'.model_pred'), sep='\t')

# Create keep files based on the results
for(i in names(pop_model_pred[,-1:-2])){
	tmp_keep<-pop_model_pred[apply(pop_model_pred[,-1:-2], 1, function(x) x[i] == max(x)),1:2]
	fwrite(tmp_keep, paste0(opt$output,'.model_pred.',i,'.keep'), sep=' ', col.names=F)
}
	
rm(targ_PCs_scaled,pop_model_pred)
gc()

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Identify individuals that are within 3SD of population specific mean for all PCs and write out scaled PCs.
###

targ_PCs<-targ_PCs[,grepl('FID|IID|PC',names(targ_PCs))]
# Read in pop_scale_for_keep
pop_scale_for_keep<-paste0(opt$output,'.',pop_keep_files$V1,'.scale')

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
	
	# Save keep file of individuals that fit the population
	fwrite(targ_PCs_scaled_i[,1:2], paste0(opt$output,'.',pop_name,'.keep'), col.names=F, sep='\t')
	
	# Write the scaled PCs
	fwrite(targ_PCs_scaled_i, paste0(opt$output,'.',pop_name,'.eigenvec'), sep='\t')
	
	rm(pop_name,pop_scale_for_keep_i,targ_PCs_scaled_i,targ_PCs_scaled_i)
	gc()
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()

