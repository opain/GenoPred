#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--target_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome target PLINK files [required]"),
make_option("--target_keep", action="store", default=NA, type='character',
		help="Path to keep file for target [optional]"),
make_option("--ref_score", action="store", default=NA, type='character',
		help="Path to reference scoring files [required]"),
make_option("--ref_freq_chr", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK .frq files [required]"),
make_option("--plink", action="store", default='plink', type='character',
		help="Path PLINK software binary [required]"),
make_option("--output", action="store", default='./Output', type='character',
		help="Path for output files [required]"),
make_option("--ref_scale", action="store", default=NA, type='character',
		help="Path reference scale file [required]"),
make_option("--pheno_name", action="store", default='./Output', type='character',
		help="Name of phenotype to be added to column names. Default is SCORE. [optional]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

if (!endsWith(opt$output_dir,'/')){
    # RM: bugfix
    opt$output_dir <- paste0(opt$output_dir, '/')
}

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# Scaled_polygenic_scorer_LDPred.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#####
# Identify shrinkage parameters
#####

ref_scale<-fread(opt$ref_scale)
param<-ref_scale$Param

#####
# Perform polygenic risk scoring
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in the target sample...')
sink()

for(param_i in param){
	if(is.na(opt$target_keep)){
		for(i in 1:22){
			system(paste0(opt$plink, ' --bfile ',opt$target_plink_chr,i,' --read-freq ',opt$ref_freq_chr,i,'.frq --score ',opt$ref_score,'.weight_',param_i,'.txt 3 4 7 sum --out ',opt$output_dir,'profiles.',param_i,'.chr',i,' --memory ',floor(opt$memory*0.9)))
		}
	} else {
		for(i in 1:22){
			system(paste0(opt$plink, ' --bfile ',opt$target_plink_chr,i,' --read-freq ',opt$ref_freq_chr,i,'.frq --keep ',opt$target_keep,' --score ',opt$ref_score,'.weight_',param_i,'.txt 3 4 7 sum --out ',opt$output_dir,'profiles.',param_i,'.chr',i,' --memory ',floor(opt$memory*0.9)))
		}
	}
}

# Add up the scores across chromosomes
profile_example<-list.files(path=opt$output_dir, pattern='*.profile')[1]
scores<-fread(paste0(opt$output_dir,profile_example))
scores<-scores[,1:2]

for(param_i in param){
	SCORE_temp<-0
	for(i in 1:22){
			profile<-fread(paste0(opt$output_dir,'profiles.',param_i,'.chr',i,'.profile'))
			SCORE_temp<-SCORE_temp+profile$SCORESUM
	}
	scores<-cbind(scores, SCORE_temp)

	names(scores)[grepl('SCORE_temp',names(scores))]<-paste0('SCORE_',param_i)
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Scale the polygenic scores based on the reference
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Scaling target polygenic scores to the reference...')
sink()

ref_scale<-fread(opt$ref_scale)
ref_scale$Param<-paste0('SCORE_',ref_scale$Param)

for(i in names(scores[,-1:-2])){
	scores[[i]]<-scores[[i]]-ref_scale$Mean[ref_scale$Param == i]
	scores[[i]]<-scores[[i]]/ref_scale$SD[ref_scale$Param == i]
	scores[[i]]<-round(scores[[i]],3)
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Write out the target sample scores
###

if(!is.na(opt$pheno_name)){
	names(scores)<-gsub('SCORE',opt$pheno_name,names(scores))
}

fwrite(scores, paste0(opt$output,'.LDPred_profiles'), sep=' ')

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Saved polygenic scores to: ',opt$output,'.LDPred_profiles.\n',sep='')
sink()

###
# Clean up temporary files
###

system(paste0('rm ',opt$output_dir,'profiles*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
