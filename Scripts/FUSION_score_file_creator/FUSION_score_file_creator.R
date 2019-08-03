#!/usr/bin/Rscript
# This script was written by Oliver Pain.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--weights", action="store", default=NA, type='character',
		help="Path for .pos file describing features [required]"),
make_option("--weights_dir", action="store", default=NA, type='character',
		help="Directory containing the weights corresponding to the features in the .pos file [required]"),
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Specify the number of cores available [required]"),
make_option("--output", action="store", default=NA, type='character',
		help="Name of output directory [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

system(paste('mkdir -p ',opt$output))

sink(file = paste(opt$output,'/FUSION_score_file_creator.log',sep=''), append = F)
cat(
'#################################################################
# FUSION_score_file_creator
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
# Convert FUSION weights into SCORE files
###################################

# Read in the .pos file
pos<-data.frame(fread(opt$weights, nThread=opt$n_cores))

sink(file = paste0(opt$output,'/FUSION_score_file_creator.log'), append = T)
cat('The .pos file contains ',dim(pos)[1],' features.\n',sep='')
sink()

# Attach weights directory to WGT values in pos file
pos$FILE<-paste0(opt$weights_dir,'/',pos$PANEL,'/',pos$PANEL,'/',sub('.*/','',pos$WGT))

# Remove .wgt.RDat from the WGT values
pos$WGT<-gsub('.wgt.RDat','',pos$WGT)
pos$WGT<-gsub('.*/','',pos$WGT)

# Create SCORE file for each set of weights using FUSION make_score.R
sink(file = paste0(opt$output,'/FUSION_score_file_creator.log'), append = T)
cat('Converting weights files into PLINK SCORE files...',sep='')
sink()

pos_chunks<-split(pos, sample(1:opt$n_cores, nrow(pos), replace=T))
rm(pos)
gc()
nrow(pos_chunks[[1]])

tmp<-foreach(chunk=1:length(pos_chunks), .combine=c) %dopar% {

	for(i in 1:nrow(pos_chunks[[chunk]])){
		# This code is the same as the FUSION make_score.R script
		load(pos_chunks[[chunk]]$FILE[i])

		best = which.min(cv.performance[2,])

		if ( names(best) == "lasso" || names(best) == "enet" ) {
		keep = wgt.matrix[,best] != 0
		} else if ( names(best) == "top1" ) {
		keep = which.max(wgt.matrix[,best]^2)
		} else { 
		keep = 1:nrow(wgt.matrix)
		}

		out<-format(cbind((snps[,c(2,5,6)]) , wgt.matrix[,best])[keep,],digits=3)
		
		fwrite(out, paste0(opt$output,'/',pos_chunks[[chunk]]$WGT[i],'.SCORE') , col.names=F , sep='\t')
		
		# Write out a snplist for each SCORE file to reduce computation time for scoring
		fwrite(out['V2'], paste0(opt$output,'/',pos_chunks[[chunk]]$WGT[i],'.snplist') , quote=F , row.names=F , col.names=F)
	
	print(i)
	}
}

sink(file = paste0(opt$output,'/FUSION_score_file_creator.log'), append = T)
cat('Done!\n',sep='')
sink()

end.time <- Sys.time()
time.taken <- end.time - start.time

sink(file = paste0(opt$output,'/FUSION_score_file_creator.log'), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
