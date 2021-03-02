#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--sumstats", action="store", default='./Output', type='character',
    help="Path for GWIZ format GWAS summary statistics [required]"),
  make_option("--gwiz", action="store", default='./Output', type='character',
    help="Path for GWIZ software [required]"),
  make_option("--output", action="store", default='./Output', type='character',
		help="Path for output files [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# gwizer.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

tmp<-sub('.*/','',opt$sumstats)
opt$sumstats_dir<-sub(paste0(tmp,'*.'),'',opt$sumstats)
opt$sumstats_file<-gsub('.*/','',opt$sumstats)
  
software_dir<-opt$gwiz
folderPath <- opt$output_dir
filename = "resampled Data"
resultsfilename = "results Data"
dir.create(paste0(opt$output_dir,filename))
dir.create(paste0(opt$output_dir,resultsfilename))
folderOut <- paste0(opt$output_dir,filename)
setwd(software_dir)

source("./modeling/GWAS_functions.0,1,2genotype (attempt 4).R")
genereteBatch(opt$sumstats_dir, folderOut, opt$sumstats_file)
study<-gsub('.csv','_Resampled.csv',opt$sumstats_file)
source("./modeling/_modelingFunctions_lda.R")
res<-batchAnalysis(paste(folderOut, "/out", sep=""), resultsfilename, study, verbose = FALSE, outerFolds = 3L, outerRep = 2L)#, classifier = "classif.logreg") 

res$aucanalysis

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('AUC =',res$aucanalysis,'\n')
sink()

unlink(paste0(opt$output_dir,filename), recursive = TRUE)
unlink(paste0(opt$output_dir,resultsfilename), recursive = TRUE)

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
