#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--keep", action="store", default=NA, type='character',
              help="File listing individuals to keep - only first column used [optional]"),
  make_option("--rel_file", action="store", default=NA, type='character',
              help="UKB relatedness file [required]"),
  make_option("--rel_thresh", action="store", default=0.044, type='numeric',
              help="Kingship threshold [optional]"),
  make_option("--seed", action="store", default=1234, type='numeric',
              help="Seed number [optional]"),
  make_option("--GreedyRelated", action="store", default=NA, type='character',
              help="Path to GreedyRelated binary [required]"),
  make_option("--output", action="store", default='./unrelated', type='character',
              help="Output file name [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat("#################################################################
# ukb_relative_remover.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at",as.character(start.time),'
Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

####
# Format UKB relatedness files for greedyRelated (modified from Joni's script)
####
library(data.table)

# Load rel_file
Related<-fread(opt$rel_file)

# Add "Pair"
Related$Pair<-1:dim(Related)[1]

# Rename "Kinship" to "Factor"
names(Related)[names(Related) == "Kinship"] <- "Factor"

# Add both IDs
Related_Long.0<-Related[,c("ID1","Pair","Factor")]
Related_Long.1<-Related[,c("ID2","Pair","Factor")]

names(Related_Long.0)[1]<-"ID"
names(Related_Long.1)[1]<-"ID"

Related_Long<-rbind(Related_Long.0, Related_Long.1)

# Sort by pairs
Related_Final<-Related_Long[order(Related_Long$Pair), ]

fwrite(Related_Final, file=paste0(opt$output,'.rel_file_temp'), sep=' ', quote=F)

######
# Now run greedyRelated
######

if(is.na(opt$keep)){
  system(paste0(opt$GreedyRelated,' -r ',opt$output,'.rel_file_temp -t ', opt$rel_thresh,' -s ',opt$seed,' > ',opt$out,'.rel_temp'))
} else {
  system(paste0(opt$GreedyRelated,' -r ',opt$output,'.rel_file_temp -keep ',opt$keep,' -t ', opt$rel_thresh,' -s ',opt$seed,' > ',opt$out,'.rel_temp'))
}

# Format output
system(paste0("cut -f 1,2 ", opt$out,'.rel_temp',' > ',opt$out,'.related'))

######
# Delete temporary file
######

system(paste0('rm ',opt$output,'.rel_file_temp'))
system(paste0('rm ',opt$output,'.rel_temp'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
