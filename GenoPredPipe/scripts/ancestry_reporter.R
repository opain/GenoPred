#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Read in samples file
config<-read.table('config.yaml', header=F)
target_list<-read.table(as.character(config$V2[config$V1 == 'target_list:']), header=T)

# Identify which super populations each target contain
ancestry_report<-NULL
target<-args[1]
output<-target_list$output[target_list$name == target]

path<-paste0(output,'/',target,'/ancestry/ancestry_all')
files<-list.files(path=path)
keep_files<-files[grepl('.keep$',files) & grepl('model_pred',files)]

info = file.info(paste0(path, '/', keep_files))
not_empty = rownames(info[info$size != 0, ])
  
ancestry_report<-rbind(ancestry_report, data.frame(name=target,
                                                   population=gsub('.keep','',gsub('.*model_pred.','',not_empty))))

write.table(ancestry_report,paste0(output,'/',target,'/ancestry/ancestry_report.txt'), col.names=T, row.names=F, quote=F)

