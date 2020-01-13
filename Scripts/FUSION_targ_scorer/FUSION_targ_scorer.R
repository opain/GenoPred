#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--targ_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome target PLINK binaries [required]"),
make_option("--targ_keep", action="store", default=NA, type='character',
		help="Path to keep file for target [optional]"),
make_option("--weights", action="store", default=NA, type='character',
		help="Path for .pos file describing features [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="RAM available in MB [required]"),
make_option("--plink", action="store", default='NA', type='character',
		help="Path to PLINK software [required]"),
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Specify the number of cores available [required]"),
make_option("--score_files", action="store", default=NA, type='character',
		help="Path to SCORE files corresponding to weights [required]"),
make_option("--ref_scale", action="store", default=NA, type='character',
		help="File containing the scale of features in the reference [required]"),
make_option("--ref_freq_chr", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK .frq files [required]"),
make_option("--pigz", action="store", default=NA, type='character',
		help="Path to pigz binary [required]"),
make_option("--output", action="store", default=NA, type='character',
		help="Name of output files [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

if(file.exists(paste0(opt$output,'.predictions.gz'))){
	print(paste0(opt$output,'.predictions.gz already exists.'))
	q()
}

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)
system(paste('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# FUSION_targ_scorer
# 4/06/2019
#################################################################

Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
registerDoMC(opt$n_cores)

# Read in the pos file
pos<-fread(opt$weights, nThread=opt$n_cores)
CHROMS<-unique(pos$CHR)[order(unique(pos$CHR))]

# Read in the ref_scale file
ref_scale<-fread(opt$ref_scale, nThread=opt$n_cores, fill=T)

sink(file = paste0(opt$output,'.log'), append = T)
cat('The .pos file contains ',dim(pos)[1],' features.\n',sep='')
sink()

# Attach weights directory to WGT values in pos file
pos$FILE<-paste0(opt$weights_dir,'/',pos$PANEL,'/',pos$PANEL,'/',sub('.*/','',pos$WGT))

# Remove .wgt.RDat from the WGT values
pos$WGT<-gsub('.wgt.RDat','',pos$WGT)
pos$WGT<-gsub('.*/','',pos$WGT)

# Create ID file to combine with predictions
if(is.na(opt$targ_keep)){
	tmp<-system(paste0(opt$plink,' --bfile ',opt$targ_plink_chr,pos$CHR[1],' --read-freq ',opt$ref_freq_chr,pos$CHR[1],'.frq --extract ',opt$score_files,'/',pos$WGT[1],'.snplist --allow-no-sex --score ',opt$score_files,'/',pos$WGT[1],'.SCORE 1 2 4 --out ',opt$output_dir,'temp.IDs --threads 1 --memory ', floor((opt$memory*0.4)/opt$n_cores)),ignore.stdout=T, ignore.stderr=T)
} else {
	tmp<-system(paste0(opt$plink,' --bfile ',opt$targ_plink_chr,pos$CHR[1],' --read-freq ',opt$ref_freq_chr,pos$CHR[1],'.frq --extract ',opt$score_files,'/',pos$WGT[1],'.snplist --keep ',opt$targ_keep,' --allow-no-sex --score ',opt$score_files,'/',pos$WGT[1],'.SCORE 1 2 4 --out ',opt$output_dir,'temp.IDs --threads 1 --memory ', floor((opt$memory*0.4)/opt$n_cores)),ignore.stdout=T, ignore.stderr=T)		
}

targ_fam<-fread(paste0(opt$output_dir,'temp.IDs.profile'), nThread=opt$n_cores)

targ_fam<-targ_fam[,1:2]
names(targ_fam)<-c('FID','IID')

write.table(targ_fam, paste0(opt$output_dir,'/TARG.IDs'), col.names=T, row.names=F, quote=F, sep=' ')
rm(targ_fam)
gc()

system(paste0('rm ',opt$output_dir,'temp.IDs.*'))

# Calculate profile scores (i.e. feature predictions)
sink(file = paste0(opt$output,'.log'), append = T)
cat('Predicting features in target sample...\n',sep='')
sink()

for(chr in CHROMS){

	sink(file = paste0(opt$output,'.log'), append = T)
	cat('Chromosome ',chr,': ',sep='')
	sink()
	
	# Subset features on the chromosome
	pos_chr<-pos[pos$CHR == chr,]
	
	sink(file = paste0(opt$output,'.log'), append = T)
	TARG_expr<-foreach(i=1:length(pos_chr$FILE), .combine=cbind) %dopar% {
	
		# Calculate feature predictions
		if(is.na(opt$targ_keep)){
			tmp<-system(paste0(opt$plink,' --bfile ',opt$targ_plink_chr,pos_chr$CHR[i],' --read-freq ',opt$ref_freq_chr,pos_chr$CHR[i],'.frq --extract ',opt$score_files,'/',pos_chr$WGT[i],'.snplist --allow-no-sex --score ',opt$score_files,'/',pos_chr$WGT[i],'.SCORE 1 2 4 --out ',opt$output_dir,pos_chr$WGT[i],' --threads 1 --memory ', floor((opt$memory*0.4)/opt$n_cores)),ignore.stdout=T, ignore.stderr=T)
		} else {
			tmp<-system(paste0(opt$plink,' --bfile ',opt$targ_plink_chr,pos_chr$CHR[i],' --read-freq ',opt$ref_freq_chr,pos_chr$CHR[i],'.frq --extract ',opt$score_files,'/',pos_chr$WGT[i],'.snplist --keep ',opt$targ_keep,' --allow-no-sex --score ',opt$score_files,'/',pos_chr$WGT[i],'.SCORE 1 2 4 --out ',opt$output_dir,pos_chr$WGT[i],' --threads 1 --memory ', floor((opt$memory*0.4)/opt$n_cores)),ignore.stdout=T, ignore.stderr=T)		
		}
		
		if(tmp != 0){ 
				# Delete temporary files
				system(paste0('rm ',opt$output_dir,pos_chr$WGT[i],'.*'),ignore.stdout=T, ignore.stderr=T)
				return(NA)
		}

		# Read in the predictions, extract SCORE column and change header.
	  feature<-fread(cmd=paste0("tr -s ' ' '\t' < ",opt$output_dir,pos_chr$WGT[i],'.profile | sed "s/^[ \t]*//" | cut -f 6'), nThread=1)
	  names(feature)<-pos_chr$WGT[i]

		# Delete temporary files
		system(paste0('rm ',opt$output_dir,pos_chr$WGT[i],'.*'),ignore.stdout=T, ignore.stderr=T)

		# Scale expression to the reference and round
		feature<-feature-ref_scale$Mean[ref_scale$WGT == pos_chr$WGT[i]]
		feature<-feature/ref_scale$SD[ref_scale$WGT == pos_chr$WGT[i]]
		feature<-round(feature,3)
		
		if(i == floor(length(pos_chr$FILE)/100*10)){cat('10% ')}
		if(i == floor(length(pos_chr$FILE)/100*20)){cat('20% ')}
		if(i == floor(length(pos_chr$FILE)/100*30)){cat('30% ')}
		if(i == floor(length(pos_chr$FILE)/100*40)){cat('40% ')}
		if(i == floor(length(pos_chr$FILE)/100*50)){cat('50% ')}
		if(i == floor(length(pos_chr$FILE)/100*60)){cat('60% ')}
		if(i == floor(length(pos_chr$FILE)/100*70)){cat('70% ')}
		if(i == floor(length(pos_chr$FILE)/100*80)){cat('80% ')}
		if(i == floor(length(pos_chr$FILE)/100*90)){cat('90% ')}
		if(i == floor(length(pos_chr$FILE)/100*100)){cat('100% ')}
		
		feature
	}

	# Remove any columns containing NA
	TARG_expr<-TARG_expr[,which(unlist(lapply(TARG_expr, function(x)!is.na(x[1])))),with=F]

	# Save TARG_expr and compress
	fwrite(TARG_expr, paste0(opt$output,'.chr',chr,'.predictions'), nThread = opt$n_cores, sep=' ')

	cat('Done!\n')
	
	sink()
	
	rm(TARG_expr)
	rm(pos_chr)
	gc()
}

# Combine per chromosome predicted expression values and insert FID and IID columns
system(paste0("paste -d ' ' $(echo ",opt$output_dir,"/TARG.IDs $(ls ",opt$output,".chr*.predictions)) | ",opt$pigz," > ",opt$output,".predictions.gz"))
system(paste0('rm ',opt$output,'.chr*.predictions'))
system(paste0('rm ',opt$output_dir,'/TARG.IDs'))

end.time <- Sys.time()
time.taken <- end.time - start.time

sink(file = paste0(opt$output,'.log'), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()



