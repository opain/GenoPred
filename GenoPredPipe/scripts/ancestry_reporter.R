#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Identify which super populations each target contain
ancestry_report<-NULL
target<-args[1]
output<-args[2]

path<-paste0(output,'/',target,'/ancestry/keep_files/model_based')
keep_files<-list.files(path=path)

# Create file list populations present in target
info = file.info(paste0(path, '/', keep_files))
not_empty = rownames(info[info$size != 0, ])
  
ancestry_report<-rbind(ancestry_report, data.frame(name=target,
                                                   population=gsub('.*\\/','',gsub('.keep','',not_empty))))

write.table(ancestry_report,paste0(output,'/',target,'/ancestry/ancestry_report.txt'), col.names=T, row.names=F, quote=F)

# Create a keep_list
keep_list<-data.frame(
    POP = gsub('.*\\/','',gsub('.keep','',not_empty)),
    file = not_empty)

write.table(keep_list, paste0(output,'/',target,'/ancestry/keep_list.txt'), col.names=T, row.names=F, quote=F)
