#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Identify which super populations each target contain
ancestry_report<-NULL
target<-args[1]
output<-args[2]

path<-paste0(output,'/',target,'/ancestry')
files<-list.files(path=path)
keep_files<-files[grepl('.keep$',files) & grepl('model_pred',files)]

info = file.info(paste0(path, '/', keep_files))
not_empty = rownames(info[info$size != 0, ])
  
ancestry_report<-rbind(ancestry_report, data.frame(name=target,
                                                   population=gsub('.keep','',gsub('.*model_pred.','',not_empty))))

write.table(ancestry_report,paste0(output,'/',target,'/ancestry/ancestry_report.txt'), col.names=T, row.names=F, quote=F)

