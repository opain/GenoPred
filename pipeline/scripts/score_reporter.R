#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Load dependencies
library(GenoUtils)
library(data.table)
source('../Scripts/functions/misc.R')

config<-args[1]

score_list <- read_param(config = config, param = 'score_list')
outdir <- read_param(config = config, param = 'outdir', return_obj = F)

if(nrow(score_list) > 0){
    score_check <- NULL
    for(name_i in score_list$name){
        log <- readLines(paste0(outdir,'/resources/data/ref/pgs_score_files/external/', name_i,'/ref-',name_i,'.log'))
        if(!any(grepl('^Skipping', log))){
            score_check <- rbind(score_check, data.frame(name=name_i, pass=T))
        } else {
            score_check <- rbind(score_check, data.frame(name=name_i, pass=F))
        }
    }
} else {
    score_check <- data.frame(
        name=NA,
        pass=NA
    )
    score_check <- score_check[complete.cases(score_check),]
}

write.table(score_check, paste0(outdir,'/resources/data/ref/pgs_score_files/external/score_report.txt'), row.names=F, quote=F)
