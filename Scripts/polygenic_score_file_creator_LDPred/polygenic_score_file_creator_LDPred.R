#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--ref_plink", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK files [required]"),
make_option("--ref_keep", action="store", default=NA, type='character',
		help="Keep file to subset individuals in reference for clumping [required]"),
make_option("--ref_pop_scale", action="store", default=NA, type='character',
		help="File containing the population code and location of the keep file [required]"),
make_option("--plink", action="store", default='plink', type='character',
		help="Path PLINK software binary [required]"),
make_option("--output", action="store", default='./Output', type='character',
		help="Path for output files [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]"),
make_option("--n_cores", action="store", default=1, type='numeric',
    help="Number of cores for parallel computing [optional]"),
make_option("--ldpred", action="store", default='ldpred', type='character',
    help="Command to call LDPred software [required]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--sumstats", action="store", default=NA, type='character',
		help="GWAS summary statistics [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

CHROMS<-1:22

if(!is.na(opt$test)){
  if(grepl('chr', opt$test)){
    single_chr_test<-T
    CHROMS<-as.numeric(gsub('chr','',opt$test))
  } else {
    single_chr_test<-F
    opt$test<-as.numeric(opt$test)
  }
}

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# polygenic_score_file_creator_LDPred.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#####
# Read in sumstats and insert p-values
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in GWAS and harmonising with reference.\n')
sink()

GWAS<-fread(cmd=paste0('zcat ',opt$sumstats))
GWAS<-GWAS[complete.cases(GWAS),]

# Extract subset if testing
if(!is.na(opt$test)){
  if(single_chr_test == F){
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Testing mode enabled. Extracted ',opt$test,' variants per chromsome.\n', sep='')
    sink()
    
    GWAS_test<-NULL
    for(i in 1:22){
      GWAS_tmp<-GWAS[GWAS$CHR == i,]
      GWAS_tmp<-GWAS_tmp[order(GWAS_tmp$BP),]
      GWAS_tmp<-GWAS_tmp[1:opt$test,]
      GWAS_test<-rbind(GWAS_test,GWAS_tmp)
    }
    
    GWAS<-GWAS_test
    GWAS<-GWAS[complete.cases(GWAS),]
    rm(GWAS_test)
    print(table(GWAS$CHR))
    
  } else {
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Testing mode enabled. Extracted chromosome ',opt$test,' variants per chromsome.\n', sep='')
    sink()
    
    GWAS<-GWAS[GWAS$CHR == CHROMS,]
    print(table(GWAS$CHR))
  }
}

# Convert OR into BETA
if(sum(names(GWAS) == 'OR') == 1){
  GWAS$BETA<-log(GWAS$OR)
}

# Rename allele frequency column
if(sum(names(GWAS) == 'FREQ') == 0){
  GWAS$FREQ<-GWAS$REF.FREQ
}

# Check GWAS_N
GWAS_N<-mean(GWAS$N)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

# Insert fake INFO score if not present
if(sum(names(GWAS) == 'INFO') == 0){
  GWAS$INFO<-1
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('INFO column not present.\n')
  sink()
}

GWAS<-GWAS[,c('CHR','BP','A1','A2','FREQ','INFO','SNP','P','BETA'),with=F]
names(GWAS)<-c('chr','pos','alt','ref','reffrq','info','rs','pval','effalt')

GWAS$chr<-paste0('chr',GWAS$chr)

fwrite(GWAS, paste0(opt$output_dir,'GWAS_sumstats.txt'), sep=' ', na = "NA", quote=F)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After harmonisation with the reference,',dim(GWAS)[1],'variants remain.\n')
sink()

if(!is.na(opt$ref_keep)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('ref_keep used to subset reference genotype data.\n')
  sink()
  
  #####
  # Create subset of ref files
  #####
  
  system(paste0(opt$plink,' --bfile ',opt$ref_plink,' --keep ',opt$ref_keep,' --make-bed --out ',opt$output_dir,'ldpred_ref'))

  opt$ref_plink_subset<-paste0(opt$output_dir,'ldpred_ref')
} else {
  opt$ref_plink_subset<-opt$ref_plink
}

if(!is.na(opt$test)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  test_start.time <- Sys.time()
  cat('Test started at',as.character(test_start.time),'\n')
  sink()
}

#####
# Coordinate the GWAS summary statistics with reference using LDPred
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Corrdinating GWAS summary statistics with reference using LDPred...')
sink()

log<-system(paste0(opt$ldpred,' coord --ssf-format STANDARD --N ',round(GWAS_N,0),' --ssf ',opt$output_dir,'GWAS_sumstats.txt --out ',opt$output_dir,'GWAS_sumstats.coord --gf ',opt$ref_plink_subset), intern=T)

writeLines(log, paste(opt$output_dir,'LDPred_coord.log',sep=''))

log_nsnp<-log[grepl("SNPs retained after filtering:", log)]
nsnp<-as.numeric(gsub(' ','', gsub('.*:','',log_nsnp)))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

#####
# Adjust the effect size estimates
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Adjusting effect sizes using LDPred...')
sink()

log2<-system(paste0(opt$ldpred,' gibbs --cf ',opt$output_dir,'GWAS_sumstats.coord --ldr ',round(nsnp/3000,0),' --ldf ', opt$ref_plink_subset,' --out ',opt$output,'.weight --N ',round(GWAS_N,0)), intern=T)

writeLines(log2, paste(opt$output_dir,'LDPred_gibbs.log',sep=''))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

if(!is.na(opt$test)){
  if(!file.exists(paste0(opt$output,'.weight_LDpred-inf.txt'))){
    q()
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  system(paste0('rm ',opt$output_dir,'*.txt'))
  system(paste0('rm ',opt$output_dir,'LDPred_gibbs.log'))
  system(paste0('rm ',opt$output_dir,'LDPred_coord.log'))
  system(paste0('rm ',opt$output_dir,'GWAS_sumstats*'))
  if(!is.na(opt$ref_keep)){
    system(paste0('rm ',opt$output_dir,'ldpred_ref*'))
  }
  q()
}

####
# Calculate mean and sd of polygenic scores at each threshold
####

# Calculate polygenic scores for reference individuals
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

score_files<-list.files(path=opt$output_dir,pattern='LDpred')
param<-gsub('.*.weight_','',score_files)
param<-gsub('.txt.*','',param)

for(i in param){
  system(paste0(opt$plink, ' --bfile ',opt$ref_plink,' --score ',opt$output,'.weight_',i,'.txt 3 4 7 sum --out ',opt$output_dir,'ref.profiles.',i,' --memory ',floor(opt$memory*0.7)))
}

fam<-fread(paste0(opt$ref_plink,'.fam'))
scores<-fam[,1:2]
names(scores)<-c('FID','IID')

for(i in param){
  SCORE_temp<-fread(paste0(opt$output_dir,'ref.profiles.',i,'.profile'))
  scores<-cbind(scores, SCORE_temp[, 6])
  names(scores)[grepl('SCORESUM',names(scores))]<-i
}

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

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Clean up temporary files
###

system(paste0('rm ',opt$output_dir,'ref.profiles.*'))
system(paste0('rm ',opt$output_dir,'GWAS_sumstats*'))
if(!is.na(opt$ref_keep)){
  system(paste0('rm ',opt$output_dir,'ldpred_ref*'))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
