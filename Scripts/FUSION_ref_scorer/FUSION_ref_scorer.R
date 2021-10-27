#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at Cardiff University under the supervision of Richard Anney.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--ref_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK binaries [required]"),
make_option("--weights", action="store", default=NA, type='character',
		help="Path for .pos file describing features [required]"),
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Specify the number of cores available [required]"),
make_option("--score_files", action="store", default=NA, type='character',
		help="Path to SCORE files corresponding to weights [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="RAM available in MB [required]"),
make_option("--plink", action="store", default='NA', type='character',
		help="Path to PLINK software [required]"),
make_option("--output", action="store", default=NA, type='character',
		help="Name of output directory [required]"),
make_option("--ref_pop_scale", action="store", default=NA, type='character',
		help="File containing the population code and location of the keep file [required]"),
make_option("--pigz", action="store", default=NA, type='character',
		help="Path to pigz binary [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

if(file.exists(paste0(opt$output,'.predictions.gz'))){
	print(paste0(opt$output,'.predictions.gz already exists.'))
	q()
}
	
opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# FUSION_ref_scorer
# 24/03/2019
#################################################################

Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
registerDoMC(opt$n_cores)

###################################
# Predict features in the reference sample
###################################

#Read in the pos file
pos<-fread(opt$weights, nThread=opt$n_cores)

sink(file = paste0(opt$output,'.log'), append = T)
cat('The .pos file contains ',dim(pos)[1],' features.\n',sep='')
sink()

# Attach weights directory to WGT values in pos file
pos$FILE<-paste0(opt$weights_dir,'/',pos$PANEL,'/',pos$PANEL,'/',sub('.*/','',pos$WGT))

# Remove .wgt.RDat from the WGT values
pos$WGT<-gsub('.wgt.RDat','',pos$WGT)
pos$WGT<-gsub('.*/','',pos$WGT)

# Read in reference fam file
ref_fam<-fread(paste0(opt$ref_plink_chr,'1.fam'), nThread=opt$n_cores)
ref_fam<-ref_fam[,1:2]
names(ref_fam)<-c('FID','IID')
	
# Calculate profile scores (i.e. feature predictions)
sink(file = paste0(opt$output,'.log'), append = T)
cat('Predicting features in reference sample...\n',sep='')
sink()

# Create column IDs to be combined to the feature predictions.
write.table(ref_fam, paste0(opt$output_dir,'REF.IDs'), col.names=T, row.names=F, quote=F, sep='\t')

error_table_all<-NULL

error_table<-foreach(i=1:length(pos$FILE), .combine=rbind) %dopar% {
	# Calculate feature predictions
	error<-system(paste0(opt$plink,' --bfile ',opt$ref_plink_chr,pos$CHR[i],' --extract ',opt$score_files,'/',pos$WGT[i],'.snplist --allow-no-sex --score ',opt$score_files,'/',pos$WGT[i],'.SCORE 1 2 4 --out ',opt$output_dir,pos$WGT[i],' --threads 1 --memory ', floor((opt$memory*0.4)/opt$n_cores)),ignore.stdout=T, ignore.stderr=T)
	# Delete temporary files and extract feature prediction column to reduce disk space
	system(paste0("awk '{print $6}' ",opt$output_dir,pos$WGT[i],'.profile | tail -n +2 > ',opt$output_dir,pos$WGT[i],'.profile_mini'),intern=T)
	system(paste0("echo ",pos$WGT[i]," | cat - ",opt$output_dir,pos$WGT[i],'.profile_mini > ',opt$output_dir,pos$WGT[i],'.profile_mini_tmp && mv ',opt$output_dir,pos$WGT[i],'.profile_mini_tmp ',opt$output_dir,pos$WGT[i],'.profile_mini'),intern=T)
	system(paste0('rm ',opt$output_dir,pos$WGT[i],'.profile'),ignore.stdout=T, ignore.stderr=T)
	system(paste0('rm ',opt$output_dir,pos$WGT[i],'.nosex'),ignore.stdout=T, ignore.stderr=T)
	system(paste0('rm ',opt$output_dir,pos$WGT[i],'.nopred'),ignore.stdout=T, ignore.stderr=T)
	print(i)
	data.frame(N=i,Error=error)
}

# Split feature predictions into list of <1000 files, then paste each list of files in batches, and then past all batches. 
system(paste0("ls -1 ", opt$output_dir,"*.profile_mini | split -l 1000 -a 4 -d - ",opt$output_dir,"profile_mini_lists"))
system(paste0("echo ", opt$output_dir,"REF.IDs | cat - ",opt$output_dir,"profile_mini_lists0000 > ",opt$output_dir,"profile_mini_lists0000_temp && mv ",opt$output_dir,"profile_mini_lists0000_temp ",opt$output_dir,"profile_mini_lists0000"))
tmp<-foreach(k=list.files(path=opt$output_dir, pattern="profile_mini_lists*"), .combine=c) %dopar% {
	system(paste0("paste $(cat ",opt$output_dir,k,") > ", opt$output_dir,"merge_",k))
}
system(paste0("paste ", opt$output_dir,"merge_profile_mini_lists* > ", opt$output,'.predictions'))
system(paste0('rm ',opt$output_dir,'REF.IDs'),ignore.stdout=T, ignore.stderr=T)

# Delete temporary files
system(paste0("rm ",opt$output_dir,'profile_mini_lists*'))
system(paste0("rm ",opt$output_dir,'merge_profile_mini_lists*'))
system(paste0("rm ",opt$output_dir,'*.profile_mini'))

# Output file containing list of features that couldn't be predicted.
error_table$Error[error_table$Error > 0]<-1
error_table<-error_table[error_table$Error > 0,]
if(dim(error_table)[1] > 0){
	for(i in 1:dim(error_table)[1]){
		k<-error_table$N[i]
		error_table$ID[i]<-pos$WGT[k]
		tmp1<-read.table(paste0(opt$output_dir,pos$WGT[k],'.log'),sep='*')
		NoValid<-sum(grepl('Error: No valid entries in --score file.',tmp1$V1))
		error_table$Reason[i]<-'Error: No valid entries in --score file.'
	}
	error_table$Error<-NULL
	error_table<-error_table[c('ID','Reason')]
	error_table_all<-rbind(error_table_all,error_table)

	write.table(error_table_all, paste0(opt$output_dir,'Prediction_failed.txt'), col.names=T, row.names=F, quote=T)
	sink(file = paste0(opt$output_dir,'.log'), append = T)
		cat(paste0(dim(error_table_all)[1],' feature/s cannot not be predicted (',opt$output_dir,'Prediction_failed.txt)\n'))
	sink()
	rm(error_table)
	rm(error_table_all)
	gc()
}

error_table<-foreach(i=1:length(pos$FILE), .combine=rbind) %dopar% {
		system(paste0('rm ',opt$output_dir,pos$WGT[i],'.log'),ignore.stdout=T, ignore.stderr=T)
}

sink(file = paste0(opt$output_dir,'.log'), append = T)
	cat('Done!\n')
sink()

###################################
# Calculate the mean and SD of feature predictions in the reference
###################################

pop_keep_files<-read.table(opt$ref_pop_scale, header=F, stringsAsFactors=F)

REF_expr<-fread(paste0(opt$output,".predictions"), nThread=opt$n_cores)

for(k in 1:dim(pop_keep_files)[1]){
	pop<-pop_keep_files$V1[k]
	keep<-fread(pop_keep_files$V2[k], header=F)
	REF_expr_keep<-REF_expr[(REF_expr$FID %in% keep$V1),]

	ref_scale<-data.frame(	WGT=names(REF_expr_keep[,-1:-2]),
													Mean=sapply(REF_expr_keep[,-1:-2], function(x) mean(x)),
													SD=sapply(REF_expr_keep[,-1:-2], function(x) sd(x)))

	fwrite(ref_scale, paste0(opt$output,'.',pop,'.scale'), sep=' ')
}

system(paste0(opt$pigz,' ',opt$output,'.predictions'))

end.time <- Sys.time()
time.taken <- end.time - start.time

sink(file = paste0(opt$output,'.log'), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
