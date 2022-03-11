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
make_option("--output", action="store", default='./Output', type='character',
		help="Path for output files [required]"),
make_option("--ref_scale", action="store", default=NA, type='character',
		help="Path reference scale file [required]"),
make_option("--covar_model", action="store", default=NA, type='character',
		help="Path to reference-derived model for covariates [optional]"),
make_option("--target_covar", action="store", default=NA, type='character',
		help="Target sample covariates data [optional]"),
make_option("--pheno_name", action="store", default='./Output', type='character',
		help="Name of phenotype to be added to column names. Default is SCORE. [optional]"),
make_option("--prsice_path", action="store", default=NA, type='character',
    help="Path to PRSice. [optional]"),
make_option("--rscript", action="store", default='Rscript', type='character',
    help="Path to Rscript [optional]"),
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
# Scaled_polygenic_scorer_dense.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#####
# Perform polygenic risk scoring
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in the target sample...')
sink()

if(is.na(opt$target_keep)){
	  system(paste0(opt$rscript,' ',opt$prsice,'/PRSice.R --prsice ',opt$prsice,'/PRSice_linux --base ',opt$ref_score,' --target ',opt$target_plink_chr,"# --thread 1 --lower 1e-8 --stat BETA --binary-target F --score sum --no-clump --no-regress --out ",opt$output,'score'))
} else {
	  system(paste0(opt$rscript,' ',opt$prsice,'/PRSice.R --prsice ',opt$prsice,'/PRSice_linux --base ',opt$ref_score,' --target ',opt$target_plink_chr,"# --keep ",opt$target_keep,' --thread 1 --lower 1e-8 --stat BETA --binary-target F --score sum --no-clump --no-regress --out ',opt$output,'score'))
}

# Read in the scores

if (file.exists(paste0(opt$output,'score.all.score'))){
    scores<-fread(paste0(opt$output,'score.all.score'))
} else {
    scores<-fread(paste0(opt$output,'score.all_score'))
}

names(scores)[-1:-2]<-paste0('SCORE_',names(scores)[-1:-2])

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

scores_scaled<-scores
for(i in names(scores)[-1:-2]){
	scores_scaled[[i]]<-scores[[i]]-ref_scale$Mean[ref_scale$pT == i]
	scores_scaled[[i]]<-scores_scaled[[i]]/ref_scale$SD[ref_scale$pT == i]
	scores_scaled[[i]]<-round(scores_scaled[[i]],3)
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Write out the target sample scores
###

if(!is.na(opt$pheno_name)){
	names(scores_scaled)<-gsub('SCORE',opt$pheno_name,names(scores_scaled))
}

fwrite(scores_scaled, paste0(opt$output,'.profiles'), sep=' ')

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Saved polygenic scores to: ',opt$output,'.profiles.\n',sep='')
sink()

###
# Clean up temporary files
###

system(paste0('rm ',opt$output,'score*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
