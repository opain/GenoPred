#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--target_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome target PLINK files [required]"),
make_option("--target_keep", action="store", default=NA, type='character',
		help="Path to keep file for target sample individuals [optional]"),
make_option("--ref_freq_chr", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK .frq files [required]"),
make_option("--ref_eigenvec", action="store", default=NA, type='character',
		help="Path to reference PLINK .eigenvec file [required]"),
make_option("--ref_score", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK .eigenvec.var file [required]"),
make_option("--plink2", action="store", default='plink', type='character',
		help="Path PLINK software binary [required]"),
make_option("--output", action="store", default='./PC_projector_output/Output', type='character',
		help="Path for output files [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]"),
make_option("--pop_data", action="store", default=NA, type='character',
		help="Path to file containing reference population data [optional]"),
make_option("--pop_model", action="store", default=NA, type='character',
		help="rds file containing population prediction model [optional]"),
make_option("--pop_model_scale", action="store", default=NA, type='character',
		help="Path to reference scaling file used when deriving the pop_model [optional]"),
make_option("--pop_scale_for_keep", action="store", default=NA, type='character',
		help="Path to population specific scale files for excluding outliers [optional]"),
make_option("--ref_scale", action="store", default=NA, type='character',
		help="Path to population specific scale file [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(ggplot2)
library(cowplot)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# scaled_ancestry_scorer.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

if(is.na(opt$pop_scale_for_keep) & is.na(opt$ref_scale)){
	sink(file = paste(opt$output,'.log',sep=''), append = F)
	cat('pop_scale_for_keep or ref_scale must be specified to scale the eigenvectors\n')
	sink()
	q()
}

#####
# Calculate PCs in target sample
#####

# Count the number of PCs in the score files.
ref_PCs<-fread(opt$ref_eigenvec)
opt$n_pcs<-dim(ref_PCs)[2]-2

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating PCs in the target sample...')
sink()

system(paste0('cut -f 2 ',opt$ref_score,' | tail -n +2 > ',opt$output_dir,'score_file.snplist'))

if(is.na(opt$target_keep)){
	for(i in 1:22){
		system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --read-freq ',opt$ref_freq_chr,i,'.frq --extract ',opt$output_dir,'score_file.snplist --score ',opt$ref_score,' header-read 2 3 --threads 1 --score-col-nums 5-',as.numeric(opt$n_pcs)+4,' --out ',opt$output_dir,'profiles.chr',i,' --memory ',floor(opt$memory*0.9)))
	}
} else {
	for(i in 1:22){
		system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --read-freq ',opt$ref_freq_chr,i,'.frq --keep ',opt$target_keep,' --extract ',opt$output_dir,'score_file.snplist --score ',opt$ref_score,' header-read 2 3 --threads 1 --score-col-nums 5-',as.numeric(opt$n_pcs)+4,' --out ',opt$output_dir,'profiles.chr',i,' --memory ',floor(opt$memory*0.9)))
	}
}
system(paste0('rm ',opt$output_dir,'score_file.snplist'))

# Add up the scores across chromosomes
scores<-fread(cmd=paste0('cut -f 1-2 ',opt$output_dir,'profiles.chr22.sscore'))
names(scores)<-c('FID','IID')

nsnp_all<-0
for(i in 1:22){
	profile<-data.frame(fread(paste0(opt$output_dir,'profiles.chr',i,'.sscore')))
	profile<-as.matrix(profile[,grepl('^PC', names(profile))])
	nsnp<-system(paste0('wc -l ',opt$ref_freq_chr,i,'.frq'), intern=T)
	nsnp<-as.numeric(unlist(strsplit(nsnp,' '))[1])
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
rm(scores,profile_all)
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

# Read in reference sample PCs
ref_PCs<-data.frame(fread(paste0(opt$ref_eigenvec)))
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

if(!is.na(opt$pop_model)){
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Estimating probability of each population...')
	sink()

	# Read in the reference scale file
	pop_model_scale<-fread(opt$pop_model_scale)

	# Scale the target PCs
	targ_PCs_scaled<-targ_PCs
	for(i in 1:dim(pop_model_scale)[1]){
		targ_PCs_scaled[[paste0('PC',i)]]<-targ_PCs[[paste0('PC',i)]]-pop_model_scale$Mean[pop_model_scale$PC == paste0('PC',i)]
		targ_PCs_scaled[[paste0('PC',i)]]<-targ_PCs_scaled[[paste0('PC',i)]]/pop_model_scale$SD[pop_model_scale$PC == paste0('PC',i)]
	}
	
	# Read in model
	pop_model<-readRDS(opt$pop_model)
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
}

###
# Identify individuals that are within 3SD of population specific mean for all PCs and write out scaled PCs.
###

if(!is.na(opt$pop_scale_for_keep)){
	targ_PCs<-targ_PCs[,grepl('FID|IID|PC',names(targ_PCs))]
	# Read in pop_scale_for_keep
	pop_scale_for_keep<-as.character(fread(opt$pop_scale_for_keep, header=F)$V1)
	for(i in 1:length(pop_scale_for_keep)){
		# Idenitfy name of population based on scale file
		pop_name<-gsub('.scale','',substr(pop_scale_for_keep[i], nchar(pop_scale_for_keep[i])-9+1, nchar(pop_scale_for_keep[i])))
		pop_scale_for_keep_i<-fread(pop_scale_for_keep[i])
		
		# Scale the target based on population scale file
		targ_PCs_scaled_i<-targ_PCs
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
}

###
# If not already done so, scale the eigenvectors
###

if(is.na(opt$pop_scale_for_keep)){
	targ_PCs<-targ_PCs[,grepl('FID|IID|PC',names(targ_PCs))]
	# Extract pop name from ref_scale file name
	pop_name<-gsub('.scale','',substr(opt$ref_scale, nchar(opt$ref_scale)-9+1, nchar(opt$ref_scale)))
	
	# Read in ref_scale
	ref_scale<-fread(opt$ref_scale)

	# Scale the target based on population scale file
	targ_PCs_scaled<-targ_PCs
	for(j in 1:dim(ref_scale)[1]){
		targ_PCs_scaled[[paste0('PC',j)]]<-targ_PCs[[paste0('PC',j)]]-ref_scale$Mean[ref_scale$PC == paste0('PC',j)]
		targ_PCs_scaled[[paste0('PC',j)]]<-targ_PCs_scaled[[paste0('PC',j)]]/ref_scale$SD[ref_scale$PC == paste0('PC',j)]
		targ_PCs_scaled[[paste0('PC',j)]]<-round(targ_PCs_scaled[[paste0('PC',j)]],3)
	}

	# Write the scaled PCs
	fwrite(targ_PCs_scaled, paste0(opt$output,'.',pop_name,'.eigenvec'), sep='\t')

	rm(pop_name,ref_keep,targ_PCs_scaled)
	gc()
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
