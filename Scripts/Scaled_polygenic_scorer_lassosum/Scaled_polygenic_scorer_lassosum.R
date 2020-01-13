#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--target_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome target PLINK files [required]"),
make_option("--target_keep", action="store", default=NA, type='character',
		help="Path to keep file for target [optional]"),
make_option("--ref_model", action="store", default=NA, type='character',
		help="Path to reference scoring files [required]"),
make_option("--output", action="store", default='./Output', type='character',
		help="Path for output files [required]"),
make_option("--ref_scale", action="store", default=NA, type='character',
		help="Path reference scale file [required]"),
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Path reference scale file [required]"),
make_option("--pheno_name", action="store", default='./Output', type='character',
		help="Name of phenotype to be added to column names. Default is SCORE. [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(lassosum)
library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# Scaled_polygenic_scorer_lassosum.R V1.0
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

# Read in the lassosum model
mod<-readRDS(opt$ref_model)

# Make fake pheno so validate runs
fam<-fread(paste0(opt$target_plink_chr,'22.fam'))
pheno<-data.frame(fam[,1:2],pheno=rnorm(dim(fam)[1]))
names(pheno)<-c('FID','IID','pheno')

is.odd <- function(x){x %% 2 != 0}
CHROMS_mat<-matrix(NA,nrow=opt$n_cores, ncol=ceiling(22/opt$n_cores)) 
CHROMS_mat[1:22]<-1:22
for(i in which(is.odd(1:dim(CHROMS_mat)[2]))){CHROMS_mat[,i]<-rev(CHROMS_mat[,i])} 
CHROMS<-as.numeric(CHROMS_mat) 
CHROMS<-CHROMS[!is.na(CHROMS)] 
print(CHROMS)

keep<-fread(opt$target_keep)
write.table(keep[,1:2], paste0(opt$output,'keep_temp.txt'), col.names=F, row.names=F, quote=F)

scores_all<-foreach(i=1:22) %dopar% {
	v2 <- validate(mod, test.bfile=paste0(opt$target_plink_chr,i), keep=paste0(opt$output,'keep_temp.txt'), pheno=pheno, plot=F)
	
	scores<-v2$results.table[!is.na(v2$results.table$order),c('FID','IID')]
	for(j in 1:length(v2$s)){
		for(k in 1:length(v2$lambda)){
			scores_tmp<-data.frame(v2$pgs[[j]][,k])
			names(scores_tmp)<-paste0('s',v2$s[j],'_lambda',v2$lambda[k])
			scores<-cbind(scores,scores_tmp)
		}
	}
	scores
}

scores_GW<-cbind(scores_all[[1]][,c('FID','IID')], Reduce('+', lapply(scores_all, "[", , -1:-2)))
rm(scores_all)
gc()

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

for(i in names(scores_GW[,-1:-2])){
	scores_GW[[i]]<-scores_GW[[i]]-ref_scale$Mean[ref_scale$Param == i]
	scores_GW[[i]]<-scores_GW[[i]]/ref_scale$SD[ref_scale$Param == i]
	scores_GW[[i]]<-round(scores_GW[[i]],3)
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Write out the target sample scores
###

if(!is.na(opt$pheno_name)){
	names(scores_GW)[-1:-2]<-paste0(opt$pheno_name,'_',names(scores_GW)[-1:-2])
}

fwrite(scores_GW, paste0(opt$output,'.lassosum_profiles'), sep=' ', na='NA')

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Saved polygenic scores to: ',opt$output,'.lassosum_profiles.\n',sep='')
sink()

system(paste0('rm ',opt$output,'keep_temp.txt'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
