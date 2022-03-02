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
make_option("--plink2", action="store", default='plink', type='character',
		help="Path PLINK v2 software binary [required]"),
make_option("--output", action="store", default='./Output', type='character',
		help="Path for output files [required]"),
make_option("--ref_scale", action="store", default=NA, type='character',
		help="Path reference scale file [required]"),
make_option("--covar_model", action="store", default=NA, type='character',
		help="Path to reference-derived model for covariates [optional]"),
make_option("--target_covar", action="store", default=NA, type='character',
		help="Target sample covariates data [optional]"),
make_option("--pheno_name", action="store", default=NA, type='character',
    help="Name of phenotype to be added to column names. Default is SCORE. [optional]"),
make_option("--n_cores", action="store", default=1, type='numeric',
    help="Number of cores to use [optional]"),
make_option("--extract", action="store", default=NA, type='character',
    help="SNP list to extract before scoring [optional]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# Scaled_polygenic_scorer.R V2.0
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

# Determine the number of scores
score_small<-fread(opt$ref_score, nrows=5)
n_scores<-ncol(score_small)-2

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Score file contains ',n_scores,' sets of coeficients.\n', sep='')
sink()

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in the target sample...')
sink()

if(n_scores > 1){
  if(is.na(opt$target_keep)){
  	for(i in 1:22){
  	  if(is.na(opt$extract)){
  		  system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --read-freq ',opt$ref_freq_chr,i,'.frq --score ',opt$ref_score,' header-read cols=fid,nallele,denom,dosagesum,scoresums --score-col-nums 3-',2+n_scores,' --out ',opt$output_dir,'profiles.chr',i,' --threads ',opt$n_cores,' --memory ',floor(opt$memory*0.9)))
  	  } else {
  	    system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --read-freq ',opt$ref_freq_chr,i,'.frq --extract ',opt$extract,' --score ',opt$ref_score,' header-read cols=fid,nallele,denom,dosagesum,scoresums --score-col-nums 3-',2+n_scores,' --out ',opt$output_dir,'profiles.chr',i,' --threads ',opt$n_cores,' --memory ',floor(opt$memory*0.9)))
  	  }
  	}
  } else {
  	for(i in 1:22){
  	  if(is.na(opt$extract)){
  	    system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --read-freq ',opt$ref_freq_chr,i,'.frq --keep ',opt$target_keep,' --score ',opt$ref_score,' header-read cols=fid,nallele,denom,dosagesum,scoresums --score-col-nums 3-',2+n_scores,' --out ',opt$output_dir,'profiles.chr',i,' --threads ',opt$n_cores,' --memory ',floor(opt$memory*0.9)))
  	  } else {
  	    system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --read-freq ',opt$ref_freq_chr,i,'.frq --extract ',opt$extract,' --keep ',opt$target_keep,' --score ',opt$ref_score,' header-read cols=fid,nallele,denom,dosagesum,scoresums --score-col-nums 3-',2+n_scores,' --out ',opt$output_dir,'profiles.chr',i,' --threads ',opt$n_cores,' --memory ',floor(opt$memory*0.9)))
  	  }
  	}
  }
} else {
  if(is.na(opt$target_keep)){
    for(i in 1:22){
      if(is.na(opt$extract)){
        system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --read-freq ',opt$ref_freq_chr,i,'.frq --score ',opt$ref_score,' header-read --out ',opt$output_dir,'profiles.chr',i,' --threads ',opt$n_cores,' --memory ',floor(opt$memory*0.9)))
      } else {
        system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --read-freq ',opt$ref_freq_chr,i,'.frq --extract ',opt$extract,' --score ',opt$ref_score,' header-read cols=fid,nallele,denom,dosagesum,scoresums --out ',opt$output_dir,'profiles.chr',i,' --threads ',opt$n_cores,' --memory ',floor(opt$memory*0.9)))
      }
    }
  } else {
    for(i in 1:22){
      if(is.na(opt$extract)){
        system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --read-freq ',opt$ref_freq_chr,i,'.frq --keep ',opt$target_keep,' --score ',opt$ref_score,' header-read cols=fid,nallele,denom,dosagesum,scoresums --out ',opt$output_dir,'profiles.chr',i,' --threads ',opt$n_cores,' --memory ',floor(opt$memory*0.9)))
      } else {
        system(paste0(opt$plink2, ' --bfile ',opt$target_plink_chr,i,' --read-freq ',opt$ref_freq_chr,i,'.frq --extract ',opt$extract,' --keep ',opt$target_keep,' --score ',opt$ref_score,' header-read cols=fid,nallele,denom,dosagesum,scoresums --out ',opt$output_dir,'profiles.chr',i,' --threads ',opt$n_cores,' --memory ',floor(opt$memory*0.9)))
      }
    }
  }  
}

# Add up the scores across chromosomes
for(i in 1:22){
  if(file.exists(paste0(opt$output_dir,'profiles.chr',i,'.sscore'))){
    scores_ids<-fread(paste0(opt$output_dir,'profiles.chr',i,'.sscore'))[,1:2, with=F]
    names(scores_ids)<-c('FID','IID')
    break
  }
}

scores<-list()
for(i in 1:22){
  if(file.exists(paste0(opt$output_dir,'profiles.chr',i,'.sscore'))){
    sscore<-fread(paste0(opt$output_dir,'profiles.chr',i,'.sscore'))
    scores[[i]]<-sscore[,grepl('SCORE_.*_SUM', names(sscore)),with=F]
    scores[[i]]<-as.matrix(scores[[i]])
  } else {
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('No scores for chromosome ',i,'. Check plink logs file for reason.\n', sep='')
    sink()
  }
}

# Remove NULL elements from list (these are inserted by R when list objects are numbered)
scores[sapply(scores, is.null)] <- NULL

scores<-Reduce(`+`, scores)
scores<-data.table(scores_ids, scores)

names(scores)<-c('FID','IID',names(score_small)[-1:-2])

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Regress out covariate effects
###

if(!is.na(opt$covar_model) == T){

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Regressing covariates based on covar_model and target_covar...')
sink()

models<-readRDS(opt$covar_model)
covar<-fread(opt$target_covar)
scores_covar<-merge(scores,covar,by=c('FID','IID'))
scores_resid<-data.frame(scores_covar[,c('FID','IID')])

for(i in names(scores[,-1:-2])){
	scores_resid[[i]]<-scores_covar[[i]]-predict(models[[i]], scores_covar)
}

scores<-scores_resid
rm(scores_resid)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()
}

###
# Scale the polygenic scores based on the reference
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Scaling target polygenic scores to the reference...')
sink()

ref_scale<-fread(opt$ref_scale)

scores_scaled<-scores
for(i in ref_scale$Param){
	scores_scaled[[i]]<-scores[[i]]-ref_scale$Mean[ref_scale$Param == i]
	scores_scaled[[i]]<-scores_scaled[[i]]/ref_scale$SD[ref_scale$Param == i]
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

system(paste0('rm ',opt$output_dir,'profiles*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
