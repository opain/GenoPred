#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--ref", action="store", default=NA, type='character',
              help="Path to bigsnpr 1000 genomes reference [required]"),
  make_option("--chr", action="store", default=1:22, type='numeric',
              help="Chromosome number [required]"),
  make_option("--pop1", action="store", default=NA, type='character',
              help="Super population code in 1KG ref [optional]"),
  make_option("--pop2", action="store", default=NA, type='character',
              help="Super population code in 1KG ref [required]"),
  make_option("--window", action="store", default=3, type='numeric',
              help="Window in cM for LD calc [required]"),
  make_option("--step", action="store", default=1, type='numeric',
              help="Step size for sliding window [required]"),
  make_option("--min_maf", action="store", default=0.01, type='numeric',
              help="Minimum MAF in both pops [required]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
              help="Number of cores available [required]"),
  make_option("--output", action="store", default="./pop_cross", type='character',
              help="Output file name [required]"),
  make_option("--memory", action="store", default=5000, type='numeric',
              help="Memory limit [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(bigsnpr)
library(bigreadr)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
  '#################################################################
# cross_pop.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

###
# Read in reference data
###

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in 1KG ref...',sep='')
sink()

# Expects format provided by bigsnpr download_1000G function
ref.bed <- bed(paste0(opt$ref,'.bed'))

if(!(file.exists(paste0(opt$ref,'.bk')))){
  snp_readBed(ref.bed$bedfile)
}
ref<-snp_attach(gsub('.bed','.rds',ref.bed$bedfile))
ref_fam_2 <- bigreadr::fread2(sub_bed(ref.bed$bedfile, ".fam2"))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n',sep='')
sink()

for(chr in opt$chr){
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('\n-------------------
Analysing chromsome ',chr,'
-------------------\n',sep='')
  sink()
  
  chr_index<-which(ref$map$chr == chr)
  
  CHR <- ref$map$chr[chr_index]
  POS <- ref$map$physical.pos[chr_index]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Calculating genetic position...',sep='')
  sink()
  
  POS2 <- snp_asGeneticPos(CHR, POS, dir ='/users/k1806347/brc_scratch/Data/Genetic_Map/CEU', ncores = opt$n_cores)
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Done!\n',sep='')
  sink()
  
  
  for(pop1 in opt$pop1){
  
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Calculating MAF and LD for ',pop1,'...',sep='')
    sink()
    
    # Identify reference individuals in pop1
    pop1_index <- which(ref_fam_2$`Super Population` == pop1)
    
    # Calculate MAF, LD, and LD-score
    pop1_maf<-bed_MAF(ref.bed, ind.row = pop1_index, ind.col = chr_index, ncores = opt$n_cores)
    pop1_corr <- snp_cor(ref$genotypes, ind.row = pop1_index, ind.col = chr_index, size = 3 / 1000, infos.pos=POS2, ncores = opt$n_cores)
    pop1_ldscore<-Matrix::colSums(pop1_corr^2)
    
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Done!\n',sep='')
    sink()
    
    for(pop2 in opt$pop2){
      
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Calculating MAF and LD for ',pop2,'...',sep='')
      sink()
      
      # Identify reference individuals in pop2
      pop2_index <- which(ref_fam_2$`Super Population` == pop2)
      
      # Calculate MAF, LD, and LD-score
      pop2_maf<-bed_MAF(ref.bed, ind.row = pop2_index, ind.col = chr_index, ncores = opt$n_cores)
      pop2_corr <- snp_cor(ref$genotypes, ind.row = pop2_index, ind.col = chr_index, size = 3 / 1000, infos.pos=POS2, ncores = opt$n_cores)
      pop2_ldscore<-Matrix::colSums(pop2_corr^2)
      
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Done!\n',sep='')
      sink()
      
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Calculating Fst between ',pop1,' and ', pop2,'...',sep='')
      sink()
      
      # Calculate the Fst
      fst_vec<-snp_fst(list(pop1_maf, pop2_maf), min_maf = 0, overall = FALSE)
      
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Done!\n',sep='')
      sink()
      
      #####
      # Calculate POPCORN C-scores (r2 covariance)
      #####
      
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Comparing LD between ',pop1,' and ', pop2,'...',sep='')
      sink()
      
      colCors = function(x, y) { 
        sqr = function(x) x*x
        x   = sweep(x, 2, colMeans(x))
        y   = sweep(y, 2, colMeans(y))
        cov = colSums(x*y)
        cor = colSums(x*y) /  sqrt(colSums(sqr(x))*colSums(sqr(y)))
        return(list(cov=cov,cor=cor))
      }
      
      # Estimate correlation using sliding window approach
      window<-opt$window
      step<-opt$step
      slider<-(2*window)+step
      
      pop1_pop2_cov_score<-NULL
      pop1_pop2_cor_score<-NULL
      j<-1
      while(j){
        print(j)
        
        start<-step*(j-1)
        stop<-step*(j-1)+slider
        
        cat('Start =',start,'\n')
        cat('Stop =',stop,'\n')
        
        block<-POS2 >= start & POS2 < stop
        
        cat('NSNP =',sum(block),'\n')
        
        tmp<-colCors(as.matrix(pop1_corr[block, block]),as.matrix(pop2_corr[block, block]))
        
        if(j == 1){
          comp_start <- 0
          comp_stop <- stop - window
        } else{
          comp_start <- start + window
          comp_stop <- stop - window
        }
        
        cat('Comp start =',comp_start,'\n')
        cat('Comp stop =',comp_stop,'\n')
        
        block_comp<-POS2[block] >= comp_start & POS2[block] < comp_stop
        
        cat('Comp NSNP =',sum(block_comp),'\n')
        
        pop1_pop2_cov_score<-c(pop1_pop2_cov_score,tmp$cov[block_comp])
        pop1_pop2_cor_score<-c(pop1_pop2_cor_score,tmp$cor[block_comp])
        
        if(comp_stop > max(POS2)){
          break()
        } else {
          j<-j+1
        }
      }
      
      sink(file = paste(opt$output,'.log',sep=''), append = T)
      cat('Done!\n',sep='')
      sink()
      
      #####
      # Tabulate data
      #####
      
      pop1_pop2_dat<-data.frame( CHR=ref$map$chromosome[chr_index],
                                 SNP=ref$map$marker.ID[chr_index],
                                 BP=ref$map$physical.pos[chr_index],
                                 A1=ref$map$allele1[chr_index],
                                 A2=ref$map$allele2[chr_index],
                                 af_1=pop1_maf$af,
                                 af_2=pop2_maf$af,
                                 Fst=fst_vec,
                                 ldscore_1=pop1_ldscore,
                                 ldscore_2=pop2_ldscore,
                                 cov_score=pop1_pop2_cov_score,
                                 cor_score=pop1_pop2_cor_score)
      
      # Remove variants with a MAF < 0.01 in either population
      pop1_pop2_dat<-pop1_pop2_dat[pop1_pop2_dat$af_1 >= 0.01 & 
                                     pop1_pop2_dat$af_1 <= 0.99 &
                                     pop1_pop2_dat$af_2 >= 0.01 & 
                                     pop1_pop2_dat$af_2 <= 0.99,]
      
      pop1_pop2_dat[,-1:-5]<-round(pop1_pop2_dat[,-1:-5],5)
      
      fwrite(pop1_pop2_dat, paste0(opt$output,'.chr',chr,'.',pop1,'.',pop2,'.txt'), col.names=T, row.names=F, quote=F, na='NA',sep='\t')
      
    }
  }
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
