#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--twas_results", action="store", default=NA, type='character',
		help="File containing TWAS results [required]"),
make_option("--ref_feature_pred", action="store", default=NA, type='character',
		help="File containing feature predictions for the reference sample [required]"),
make_option("--output", action="store", default=NA, type='character',
		help="Name of output files [required]"),
make_option("--clump_thresh", action="store", default=0.9, type='numeric',
		help="r2 threshold for clumping [optional]"),
make_option("--cor_window", action="store", default=5e6, type='numeric',
		help="Window for clumping [optional]"),
make_option("--pTs", action="store", default='1,5e-1,1e-1,5e-2,1e-2,1e-3,1e-4,1e-5,1e-6', type='character',
		help="Window for deriving pruning blocks in bases[optional]"),
make_option("--clump_mhc", action="store", default=T, type='logical',
		help="Retain only the most significant gene within the MHC region [optional]"),
make_option("--ref_keep", action="store", default=NA, type='character',
		help="Keep file for reference individuals [optional]"),
make_option("--panel", action="store", default=NA, type='character',
		help="Panel from TWAS [optional]"),
make_option("--r2_weighted", action="store", default=F, type='logical',
		help="Set to T if gene expression should be weighted by R2 of predicted expression [optional]"),
make_option("--ref_scale", action="store", default=NA, type='character',
		help="Path to file for scaling feature predictions [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

opt$pTs<- as.numeric(unlist(strsplit(opt$pTs,',')))

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)

if(file.exists(paste(opt$output,'.profiles',sep=''))){
	cat('A file named',paste(opt$output,'.profiles',sep=''),'already exists.\n')
	q()
}

system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# ref_functionally_informed_risk_scorer.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

print(opt)
cat('Analysis started at',as.character(start.time),'\n')

sink()
suppressMessages(library(data.table))

###
# Read in the input files and find intersect
###

# Read in the TWAS results
if(substr(opt$twas_results, nchar(opt$twas_results)-2, nchar(opt$twas_results)) == '.gz'){
	TWAS<-fread(paste('zcat ',opt$twas_results,sep=''))
} else {
	TWAS<-fread(opt$twas_results)
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('twas_results contains', dim(TWAS)[1],'rows.\n')

if(opt$r2_weighted == T){
###
# Weight TWAS.Z by MODELCV.R2
###

TWAS$TWAS.Z<-TWAS$TWAS.Z*TWAS$MODELCV.R2
# Recalculate TWAS.P
TWAS$TWAS.P<-2*pnorm(-abs(TWAS$TWAS.Z))

}

# Extract the weights from panel
if(!is.na(opt$panel)){
	TWAS<-TWAS[TWAS$PANEL == opt$panel,]
}

# Convert TWAS FILE column to match the gene expression column names in expression_ref
TWAS$FILE<-sub(".*/", "", TWAS$FILE)
TWAS$FILE<-sub(".wgt.RDat", "", TWAS$FILE)
TWAS$FILE<-gsub(":", ".", TWAS$FILE)
TWAS$FILE<-gsub("-", ".", TWAS$FILE)

# Remove rows with containing NAs or duplicate FILE ID.
TWAS<-TWAS[!is.na(TWAS$TWAS.Z),]
TWAS<-TWAS[!duplicated(TWAS$FILE),]

cat('twas_results contains', dim(TWAS)[1],'unique features with non-missing TWAS.Z values.\n')
if(dim(TWAS)[1] == 0){
	cat('Quitting.\n')
	q()
}
	
cat('The minimum TWAS.P value is ', min(TWAS$TWAS.P),'.\n',sep='')

# Sort the TWAS results based on CHR and P0
TWAS<-TWAS[order(TWAS$CHR,TWAS$P0),]

# Add the 0.5 Mb window to P0 and P1 values to account for window size when deriving weights.
TWAS$P0<-TWAS$P0-5e5
TWAS$P0[TWAS$P0 < 0]<-0
TWAS$P1<-TWAS$P1+5e5

sink()

# Read in the reference feature predictions
if(substr(opt$ref_feature_pred, nchar(opt$ref_feature_pred)-16, nchar(opt$ref_feature_pred)) == '.predictions_list'){
	RefGene_list<-fread(opt$ref_feature_pred, header=F)$V1
	for(i in 1:length(RefGene_list)){
		if(substr(RefGene_list[i], nchar(RefGene_list[i])-2, nchar(RefGene_list[i])) == '.gz'){
			RefGene_tmp<-fread(cmd=paste('zcat ',RefGene_list[i],sep=''))
		} else {
			RefGene_tmp<-fread(RefGene_list[i])
		}
		
		if(!is.na(opt$ref_keep)){
		# Extract individuals in the keep file
			ref_keep<-fread(opt$ref_keep, header=F)
			RefGene_tmp<-RefGene_tmp[(RefGene_tmp$FID %in% ref_keep$V1),]			
		}

		if(i == 1){
			RefGene<-RefGene_tmp
		} else {
			RefGene<-merge(RefGene, RefGene_tmp, by=c('FID','IID'))
		}
		print(i)
	}
} else {
	if(substr(opt$ref_feature_pred, nchar(opt$ref_feature_pred)-2, nchar(opt$ref_feature_pred)) == '.gz'){
		RefGene<-fread(paste('zcat ',opt$ref_feature_pred,sep=''))
	} else {
		RefGene<-fread(opt$ref_feature_pred)
	}
	
	if(!is.na(opt$ref_keep)){
	# Extract individuals in the keep file
		ref_keep<-fread(opt$ref_keep, header=F)
		RefGene<-RefGene[(RefGene$FID %in% ref_keep$V1),]
	}
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat(dim(RefGene)[1], 'individuals remain after restricting to ref_keep.\n')
sink()

# Remove features that contain NA
incomp_feat<-names(RefGene)[-1:-2][colSums(is.na(RefGene[,-1:-2])) != 0]
if(length(incomp_feat) != 0){
	RefGene<-RefGene[,!(names(RefGene) %in% incomp_feat),with=F]
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat(dim(RefGene)[2]-2, 'complete features are present in ref_feature_pred.\n')
sink()

# Convert '-' and ':' to '.' to match TWAS$FILE
names(RefGene)<-gsub('-','.',names(RefGene))
names(RefGene)<-gsub(':','.',names(RefGene))

# Extract intersecting genes between the reference feature predictions and TWAS results
intersecting_genes<-intersect(TWAS$FILE, names(RefGene))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat(length(intersecting_genes), 'features are present in both twas_results and ref_feature_pred.\n')
sink()

TWAS<-TWAS[(TWAS$FILE %in% intersecting_genes),]
RefGene<-RefGene[,c('FID','IID',intersecting_genes), with=FALSE]

if(opt$clump_mhc == T){
###
# Clump the MHC region to only contain the most significant gene in that region
###
	TWAS_notMHC<-TWAS[!(TWAS$CHR == 6 & TWAS$P0 > 26e6 & TWAS$P1 < 34e6),]
	TWAS_MHC<-TWAS[TWAS$CHR == 6 & TWAS$P0 > 26e6 & TWAS$P1 < 34e6,]
	if(dim(TWAS_MHC)[1] > 0){
		TWAS_MHC_retain<-TWAS_MHC[TWAS_MHC$TWAS.P == min(TWAS_MHC$TWAS.P),]
		TWAS<-rbind(TWAS_notMHC,TWAS_MHC_retain)
		RefGene<-RefGene[,c('FID','IID',TWAS$FILE),with=F]
		sink(file = paste(opt$output,'.log',sep=''), append = T)
			cat(dim(TWAS)[1], 'after clumping the MHC region to contain only the top gene.\n')
		sink()
	} else {
		sink(file = paste(opt$output,'.log',sep=''), append = T)
			cat('No features present in MHC region.\n')
		sink()
  }
}

###
# Create a table showing the number of genes in each pT
###
# The number of features after clumping will be filled in as each pT is processed.

NGenes_table<-NULL
for(i in 1:length(opt$pTs)){
	NGenes<-data.frame(	pT=opt$pTs[i],
						NGenes=sum(TWAS$TWAS.P < opt$pTs[i]),
						NGenes_post_clump=NA)
	NGenes_table<-rbind(NGenes_table,NGenes)
}

###
# Determine gene blocks.
###
TWAS$Block<-NA
for(i in 1:dim(TWAS)[1]){
if(i == 1){
		TWAS$Block[i]<-1
		} else {
		if(i > 1 & TWAS$CHR[i] == TWAS$CHR[i-1] & TWAS$P1[i] > (TWAS$P0[i-1] - opt$cor_window) & TWAS$P0[i] < (TWAS$P1[i-1] + opt$cor_window)){
				TWAS$Block[i]<-TWAS$Block[i-1]
				}
		if(!(i > 1 & TWAS$CHR[i] == TWAS$CHR[i-1] & TWAS$P1[i] > (TWAS$P0[i-1] - opt$cor_window) & TWAS$P0[i] < (TWAS$P1[i-1] + opt$cor_window))){
				TWAS$Block[i]<-TWAS$Block[i-1]+1
				}
		}
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('The features were be separated into',length(unique(TWAS$Block)),'blocks.\n')
sink()

###
# Clump genes based on their correlation in the FUSION 1KG reference
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Clumping features based on correlation in ref_feature_pred...')
sink()

TWAS_clumped<-NULL
for(j in unique(TWAS$Block)){
	if(sum(TWAS$Block == j) == 1){
		TWAS_clumped<-rbind(TWAS_clumped,TWAS[TWAS$Block == j,])
	} else {
		TWAS_Block<-TWAS[TWAS$Block == j,]
		TWAS_Block<-TWAS_Block[order(TWAS_Block$TWAS.P),]
		cor_block<-cor(as.matrix(RefGene[,TWAS_Block$FILE,with=F]), method='pearson')

		i<-1
		while(i){
		  # Subset rows within range of variable i
		  TWAS_Block_i<-TWAS_Block[TWAS_Block$P0 < (TWAS_Block$P1[i] + opt$cor_window) & TWAS_Block$P1 > (TWAS_Block$P0[i] - opt$cor_window),]
		  cor_block_i<-cor_block[(dimnames(cor_block)[[1]] %in% TWAS_Block_i$FILE),(dimnames(cor_block)[[2]] %in% TWAS_Block_i$FILE)]
		  
		  # If no nearby features (after previous pruning) skip
		  if(dim(TWAS_Block_i)[1] == 1){
			i<-i+1
		  	next()
		  }
		  
		  # Skip if no variables in range have correlation greater than threshold
		  if(!(max(abs(cor_block_i[-1,1])) > sqrt(opt$clump_thresh))){
		    i<-i+1
				if(i > dim(TWAS_Block)[1]){
					break()
				}
		    next()
		  }
		  
		  # Create list of variables that correlate to highly with variable i
		  exclude<-dimnames(cor_block_i)[[1]][-1][abs(cor_block_i[-1,1]) > sqrt(opt$clump_thresh)]
		  
		  # Remove these correlated variables
		  TWAS_Block<-TWAS_Block[!(TWAS_Block$FILE %in% exclude),]
		  cor_block<-cor_block[!(dimnames(cor_block)[[1]] %in% exclude),!(dimnames(cor_block)[[2]] %in% exclude)]

			i<-i+1

			# Break if we have gone through all variables.
			if(i > dim(TWAS_Block)[1]){
				break()
			}
		}	
	TWAS_clumped<-rbind(TWAS_clumped,TWAS_Block)
	}
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

# Save list of clumped TWAS results with TWAS.Z and TWAS.P values
fwrite(TWAS_clumped[,c('FILE','TWAS.Z','TWAS.P'), with=F], paste0(opt$output,'.score'), sep=' ')

# Update NGenes_table after clumping
for(i in 1:length(opt$pTs)){
	NGenes_table$NGenes_post_clump[i]<-sum(TWAS_clumped$TWAS.P <= opt$pTs[i])
}

write.table(NGenes_table, paste(opt$output,'.NFeat',sep=''), row.names=F, quote=F, sep='\t')

###
# Calculate risk scores in the reference sample for scaling the target sample scores
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating functional risk scores in ref_feature_pred...')
sink()

# Extract features after clumping
RefGene<-data.frame(RefGene[,c('FID','IID',TWAS_clumped$FILE),with=F])

# Scale the gene expression in the reference sample (should this be done using ancestry relevent scale files?)
if(substr(opt$ref_feature_pred, nchar(opt$ref_feature_pred)-16, nchar(opt$ref_feature_pred)) == '.predictions_list'){
	ref_scale_list<-fread(opt$ref_scale, header=F)$V1
	ref_scale<-NULL
	for(i in 1:length(ref_scale_list)){
		ref_scale<-rbind(ref_scale,fread(ref_scale_list[i], fill=T))
	}
} else {
	ref_scale<-fread(opt$ref_scale, fill=T)
	ref_scale$WGT<-gsub('-','.',ref_scale$WGT)
	ref_scale$WGT<-gsub(':','.',ref_scale$WGT)
}

ref_scale<-ref_scale[match(names(RefGene)[-1:-2], ref_scale$WGT),]
RefGene_ID<-RefGene[,1:2]
RefGene_noID<-t(t(RefGene[,-1:-2]) - ref_scale$Mean)
RefGene_noID<-t(t(RefGene_noID) / ref_scale$SD)
RefGene_noID<-round(RefGene_noID,3)
RefGene<-cbind(RefGene_ID,RefGene_noID)

# Weight the gene expression in each individuals by TWAS.Z
RefGene_ID<-RefGene[,1:2]
TWAS_clumped<-TWAS_clumped[match(names(RefGene)[-1:-2], TWAS_clumped$FILE),]
RefGene_noID<-t(t(RefGene[,-1:-2]) * TWAS_clumped$TWAS.Z)
RefGene<-cbind(RefGene_ID,RefGene_noID)

# For each pT calculate GeRS in reference sample using pruned gene list
GeneX_Risk_reference<-data.frame(	FID=RefGene$FID,
																	IID=RefGene$IID)

for(i in 1:sum(NGenes_table$pT > min(TWAS_clumped$TWAS.P))){
	tmp<-rowSums(RefGene[which(names(RefGene) %in% TWAS_clumped$FILE[TWAS_clumped$TWAS.P <= NGenes_table$pT[i]])])
	GeneX_Risk_reference<-cbind(GeneX_Risk_reference,tmp)
	names(GeneX_Risk_reference)[2+i]<-paste0('SCORE_',NGenes_table$pT[i])
}

# Calculate the mean and sd for each SCORE
SCORE_ref_scale<-data.frame(pT=names(GeneX_Risk_reference[-1:-2]),
							Mean=sapply(GeneX_Risk_reference[,-1:-2], function(x) mean(x)),
							SD=sapply(GeneX_Risk_reference[,-1:-2], function(x) sd(x)))

write.table(SCORE_ref_scale, paste(opt$output,'.scale',sep=''), row.names=F, quote=F, sep='\t')

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()