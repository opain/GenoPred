#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')

config<-args[1]

gwas_groups <- read_param(config = config, param = 'gwas_groups')

if(!is.null(gwas_groups)){
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  if(nrow(gwas_groups) > 0){
      leopard_check <- NULL
      for(name_i in gwas_groups$name){
          if(file.exists(paste0(outdir,'/reference/pgs_score_files/leopard/', name_i,'/ref-',name_i,'.weights.rds'))){
            leopard_check <- rbind(leopard_check, data.frame(name=name_i, pass=T))
          } else {
            leopard_check <- rbind(leopard_check, data.frame(name=name_i, pass=F))
          }
      }
  } else {
    leopard_check <- data.frame(
          name=NA,
          pass=NA
      )
    leopard_check <- leopard_check[complete.cases(leopard_check),]
  }

  write.table(leopard_check, paste0(outdir,'/reference/pgs_score_files/leopard/leopard_report.txt'), row.names=F, quote=F)
}
