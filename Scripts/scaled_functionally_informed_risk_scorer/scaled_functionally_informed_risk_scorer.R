#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--targ_feature_pred", action="store", default=NA, type='character',
		help="Path to predicted expression file [required]"),
make_option("--target_keep", action="store", default=NA, type='character',
    help="Path to keep file for target [optional]"),
make_option("--ref_score", action="store", default=NA, type='character',
		help="Path to reference scoring files [required]"),
make_option("--output", action="store", default='./Output', type='character',
		help="Path for output files [required]"),
make_option("--ref_scale", action="store", default=NA, type='character',
		help="Path reference scale file [required]"),
make_option("--pheno_name", action="store", default=NA, type='character',
		help="Name of phenotype to be added to column names. Default is SCORE.[optional]"),
make_option("--pigz", action="store", default=NA, type='character',
		help="Path to pigz binary [required]"),
make_option("--batch_size", action="store", default=5000, type='numeric',
		help="Number of individuals to be processed in each batch [optional]"),
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Specify the number of cores available [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# scaled_funcionally_informed_risk_scorer.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

# Read in the target feature predictions names
if(substr(opt$targ_feature_pred, nchar(opt$targ_feature_pred)-2, nchar(opt$targ_feature_pred)) == '.gz'){
	originally_compressed<-T
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Decompressing targ_feature_pred...')
	sink()
	system(paste0(opt$pigz,' -dc -p ',opt$n_cores,' ',opt$targ_feature_pred,' > ',opt$output,'.target_feature_pred'))
	opt$targ_feature_pred<-paste0(opt$output,'.target_feature_pred')
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Done!\n')
	sink()
} else {
originally_compressed<-F
}

TargGene<-fread(cmd=paste0('head -n 1 ',opt$targ_feature_pred), nThread=opt$n_cores)
n_indiv<-system(paste0('wc -l ', opt$targ_feature_pred), intern=T)
n_indiv<-as.numeric(gsub(' .*','',n_indiv))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat(dim(TargGene)[2]-2, 'features are present in targ_feature_pred.\n')
cat(n_indiv-1, 'individulas are present in targ_feature_pred.\n')
sink()

if(!is.na(opt$target_keep)){
  keep<-fread(opt$target_keep)
  names(keep)[1:2]<-c('FID','IID')
  write.table(keep[,1:2], paste0(opt$output,'.keep'), col.names=T, row.names=F, quote=F)
  
  system(paste0("awk -F' ' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' ",opt$output,'.keep ' ,opt$targ_feature_pred,' > ', opt$output,'.target_feature_pred.subset'))
  
  if(originally_compressed == T){
    system(paste0('rm ',opt$targ_feature_pred))
  }
  
  opt$targ_feature_pred<-paste0(opt$output,'.target_feature_pred.subset')

  n_indiv<-system(paste0('wc -l ', opt$targ_feature_pred), intern=T)
  n_indiv<-as.numeric(gsub(' .*','',n_indiv))
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat(n_indiv-1, 'individuals remain after applying keep file.\n')
  sink()
  
  system(paste0('rm ',opt$output,'.keep'))
}

# Convert '-' and ':' to '.' to match TWAS$FILE
names(TargGene)<-gsub('-','.',names(TargGene))
names(TargGene)<-gsub(':','.',names(TargGene))

# Read in the score file
TWAS_clumped<-fread(opt$ref_score, nThread=opt$n_cores)
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat(length(TWAS_clumped$FILE), 'features are present in ref_score.\n')
if(sum(names(TWAS_clumped) == 'COLOC.PP4') == 1){
  cat('COLOC PP4 estimates will be used to filter features.\n')
  opt$coloc<-T
} else {
  opt$coloc<-F
}
sink()

# Extract intersecting genes between the target feature predictions and the score file
intersecting_genes<-intersect(TWAS_clumped$FILE, names(TargGene))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat(length(intersecting_genes), 'features are present in both ref_score and targ_feature_pred.\n')
sink()

# Extract intersecting features from TWAS_clumped and sort to be the same as TargGene
TWAS_clumped<-TWAS_clumped[(TWAS_clumped$FILE %in% intersecting_genes),]
TWAS_clumped<-setDT(TWAS_clumped, key = "FILE")[names(TargGene[,intersecting_genes, with=FALSE])]

# Read in the scale file
ref_scale<-fread(opt$ref_scale)
ref_scale$pT_num<-as.numeric(gsub('SCORE_','',ref_scale$pT))
ref_scale<-ref_scale[order(ref_scale$pT_num),]

# Calculate risk scores in batches to reduce memory requirements
batch_list<-split(2:n_indiv, ceiling(seq_along(2:n_indiv)/opt$batch_size))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Data will be split into',length(batch_list),'batches.\n')
sink()

# For each pT calculate GeRS in target sample
GeneX_Risk<-foreach(batch=1:length(batch_list), .combine=rbind) %dopar% {

	# Read in the target feature predictions names
	TargGene_batch<-fread(cmd=paste0("sed -n '",min(batch_list[[batch]]),",",max(batch_list[[batch]]),"p;",max(batch_list[[batch]])+1,"q' ",opt$targ_feature_pred), nThread=1)
	names(TargGene_batch)<-names(TargGene)

	TargGene_batch<-cbind(TargGene_batch[,1:2], t(t(TargGene_batch[,intersecting_genes, with=FALSE]) * TWAS_clumped$TWAS.Z))

	GeneX_Risk_tmp<-data.table(	FID=TargGene_batch$FID,
															IID=TargGene_batch$IID)

	if(opt$coloc == T){
	  for(i in 1:sum(ref_scale$pT_num < max(TWAS_clumped$COLOC.PP4))){
	    tmp<-rowSums(TargGene_batch[,which(names(TargGene_batch) %in% TWAS_clumped$FILE[TWAS_clumped$COLOC.PP4 >= ref_scale$pT_num[i]]), with=F])
	    tmp<-tmp-ref_scale$Mean[i]
	    tmp<-tmp/ref_scale$SD[i]
	    GeneX_Risk_tmp<-cbind(GeneX_Risk_tmp,round(tmp,4))
	    names(GeneX_Risk_tmp)[2+i]<-paste0('SCORE_',ref_scale$pT_num[i])
	  }
	} else {
	  for(i in 1:sum(ref_scale$pT_num > min(TWAS_clumped$TWAS.P))){
	    tmp<-rowSums(TargGene_batch[,which(names(TargGene_batch) %in% TWAS_clumped$FILE[TWAS_clumped$TWAS.P <= ref_scale$pT_num[i]]), with=F])
	    tmp<-tmp-ref_scale$Mean[i]
	    tmp<-tmp/ref_scale$SD[i]
	    GeneX_Risk_tmp<-cbind(GeneX_Risk_tmp,round(tmp,4))
	    names(GeneX_Risk_tmp)[2+i]<-paste0('SCORE_',ref_scale$pT_num[i])
	  }
	}

	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Processed batch ',batch,'.\n', sep='')
	sink()

	print(batch)
	
	GeneX_Risk_tmp
	
}

if(dim(GeneX_Risk)[1] != (n_indiv-1)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Error during scoring.\n')
  sink()
  q()
}

names(GeneX_Risk)<-gsub('SCORE_',paste0(opt$pheno_name,'_'),names(GeneX_Risk))

fwrite(GeneX_Risk, paste0(opt$output,'.fiprofile'), sep=' ', nThread=opt$n_cores)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Finished!\n')
sink()

if(originally_compressed == T){
	system(paste0('rm ',opt$targ_feature_pred))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
