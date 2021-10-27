#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--ref_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK files [required]"),
make_option("--ref_pop_scale", action="store", default=NA, type='character',
		help="File containing the population code and location of the keep file [required]"),
make_option("--plink2", action="store", default='plink', type='character',
    help="Path PLINKv2 software binary [required]"),
make_option("--output", action="store", default='./Output', type='character',
		help="Path for output files [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]"),
make_option("--score", action="store", default=NA, type='character',
		help="Score file with format SNP A1 and then one or more effect sizes [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

CHROMS<-1:22

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# external_score_processor.R V2.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

####
# Calculate mean and sd of scores
####

score<-fread(opt$score, nrows=5)
n_score<-ncol(score)-2

# Calculate polygenic scores for reference individuals
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

if(n_score == 1){
  for(i in CHROMS){
    system(paste0(opt$plink2, ' --bfile ',opt$ref_plink_chr,i,' --score ',opt$score,' header-read --out ',opt$output,'.profiles.chr',i,' --threads 1 --memory ',floor(opt$memory*0.7)))
  }
} else {
  for(i in CHROMS){
    system(paste0(opt$plink2, ' --bfile ',opt$ref_plink_chr,i,' --score ',opt$score,' header-read --score-col-nums 3-',2+n_score,' --out ',opt$output,'.profiles.chr',i,' --threads 1 --memory ',floor(opt$memory*0.7)))
  }
}

# Add up the scores across chromosomes
fam<-fread(paste0(opt$ref_plink_chr,'22.fam'))

scores<-list()
for(i in as.character(CHROMS)){
  sscore<-fread(paste0(opt$output,'.profiles.chr',i,'.sscore'))
  scores[[i]]<-sscore[,grepl('SCORE_', names(sscore)),with=F]
  scores[[i]]<-as.matrix(scores[[i]]*sscore$NMISS_ALLELE_CT)
}

scores<-Reduce(`+`, scores)
scores<-data.table(FID=fam$V1,
                   IID=fam$V2,
                   scores)

names(scores)<-c('FID','IID',names(score)[-1:-2])

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

# Calculate the mean and sd of scores for each population specified in pop_scale
pop_keep_files<-read.table(opt$ref_pop_scale, header=F, stringsAsFactors=F)

for(k in 1:dim(pop_keep_files)[1]){
	pop<-pop_keep_files$V1[k]
	keep<-fread(pop_keep_files$V2[k], header=F)
	scores_keep<-scores[(scores$FID %in% keep$V1),]

	ref_scale<-data.frame(	Param=names(scores_keep[,-1:-2]),
													Mean=round(sapply(scores_keep[,-1:-2], function(x) mean(x)),3),
													SD=round(sapply(scores_keep[,-1:-2], function(x) sd(x)),3))

	fwrite(ref_scale, paste0(opt$output,'.',pop,'.scale'), sep=' ')
}

###
# Clean up temporary files
###

system(paste0('rm ',opt$output,'.profiles.*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
