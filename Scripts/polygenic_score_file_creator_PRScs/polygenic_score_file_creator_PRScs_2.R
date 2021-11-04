#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--ref_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK files [required]"),
make_option("--ref_keep", action="store", default=NA, type='character',
		help="Keep file to subset individuals in reference for clumping [required]"),
make_option("--ref_pop_scale", action="store", default=NA, type='character',
		help="File containing the population code and location of the keep file [required]"),
make_option("--plink", action="store", default='plink', type='character',
		help="Path PLINK software binary [required]"),
make_option("--output", action="store", default='NA', type='character',
		help="Path for output files [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]"),
make_option("--sumstats", action="store", default=NA, type='character',
		help="GWAS summary statistics [optional]"),
make_option("--PRScs_path", action="store", default=NA, type='character',
		help="Path to PRScs executable [required]"),
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Number of cores for parallel computing [optional]"),
make_option("--PRScs_ref_path", action="store", default=T, type='character',
		help="Path to PRScs reference [required]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--phi_param", action="store", default='auto', type='character',
		help="Path to PRScs reference [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

phi_param<-unlist(strsplit(opt$phi_param,','))

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

if(any(grepl('auto',phi_param))){
	phi_param<-c(sprintf("%1.00e", as.numeric(phi_param[!grepl('auto',phi_param)])),'auto')
} else {
	phi_param<-sprintf("%1.00e", as.numeric(phi_param))
}

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# polygenic_score_file_creator_PRScs.R V1.0
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

GWAS<-fread(cmd=paste0('zcat ',opt$sumstats), nThread=1)
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

GWAS_N<-mean(GWAS$N)

if(('BETA' %in% names(GWAS))){
  GWAS<-GWAS[,c('SNP','A1','A2','BETA','P')]
} else {
  GWAS<-GWAS[,c('SNP','A1','A2','OR','P')]
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

fwrite(GWAS, paste0(opt$output_dir,'GWAS_sumstats_temp.txt'), sep=' ')

rm(GWAS)
gc()

if(!is.na(opt$test)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  test_start.time <- Sys.time()
  cat('Test started at',as.character(test_start.time),'\n')
  sink()
}

#####
# Process sumstats using PRScs
#####

# Make a data.frame listing chromosome and phi combinations
jobs<-NULL
for(i in CHROMS){
  jobs<-rbind(jobs,data.frame(CHR=i,
                       phi=phi_param))
}

write.table(jobs, paste0(opt$output_dir,'job_list'), col.names=F, row.names=F, quote=F)

# Write batch job
writeLines(paste0("#!/bin/sh

#SBATCH -p shared,brc
#SBATCH --mem 5G
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH -t 24:00:00
#SBATCH -J PRScs

export MKL_NUM_THREADS=$SLURM_CPUS_ON_NODE
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_ON_NODE
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

echo $SLURM_CPUS_ON_NODE
chr=$(awk -v var=$SLURM_ARRAY_TASK_ID 'NR == var {print $1}' ", opt$output_dir,"job_list)
phi=$(awk -v var=$SLURM_ARRAY_TASK_ID 'NR == var {print $2}' ", opt$output_dir,"job_list)

echo ${chr}
echo ${phi}

if [ ${phi} == \"auto\" ];then

",opt$PRScs_path," --ref_dir=",opt$PRScs_ref_path," --bim_prefix=",opt$ref_plink_chr,"${chr}  --sst_file=",opt$output_dir,"/GWAS_sumstats_temp.txt --n_gwas=",round(GWAS_N,0)," --out_dir=",opt$output," --chrom=${chr}

else

",opt$PRScs_path," --ref_dir=",opt$PRScs_ref_path," --bim_prefix=",opt$ref_plink_chr,"${chr} --phi=${phi} --sst_file=",opt$output_dir,"/GWAS_sumstats_temp.txt --n_gwas=",round(GWAS_N,0)," --out_dir=",opt$output," --chrom=${chr}

fi

"), paste0(opt$output_dir,'batch.sh'))

# Run batch job
jobID<-system(paste0("sbatch --array ",1,"-",nrow(jobs),"%",opt$n_cores," ", opt$output_dir,'batch.sh'),intern=T) 
jobID<-gsub('.* ','', jobID)

# Check whether finished
Sys.sleep(30)
while(i){
  system(paste0('sacct -j ',jobID,' > ',opt$output_dir,'sacct_log.txt'))
  sacct_log<-fread(paste0(opt$output_dir,'sacct_log.txt'), fill=T)
  sacct_log<-sacct_log[sacct_log$JobName == 'PRScs',]
  
  print(sacct_log)
  
  if(sum(sacct_log$State == 'FAILED') > 0){
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Job failed.\n')
    sink()
    q()
  }
  
  if(sum(sacct_log$State != 'COMPLETED') == 0){
    break
  } else {
    Sys.sleep(60)
  }
}

# Once finished, delete temporary GWAS_stats_temp file
system(paste0('rm ',opt$output_dir,'GWAS_sumstats_temp.txt'))

if(!is.na(opt$test)){
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  system(paste0('rm ',opt$output,'*.txt'))
  q()
}

####
# Calculate mean and sd of polygenic scores
####

# Calculate polygenic scores for reference individuals
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

for(phi_i in phi_param){
	for(i in CHROMS){
		system(paste0(opt$plink, ' --bfile ',opt$ref_plink_chr,i,' --score ',opt$output,'_pst_eff_a1_b0.5_phi',phi_i,'_chr',i,'.txt 2 4 6 sum --out ',opt$output,'.',phi_i,'.profiles.chr',i,' --memory ',floor(opt$memory*0.7)))
	}
}

# Add up the scores across chromosomes
fam<-fread(paste0(opt$ref_plink_chr,'22.fam'))
scores<-fam[,1:2]
names(scores)<-c('FID','IID')

for(phi_i in phi_param){

	SCORE_temp<-0
	for(i in CHROMS){
			profile<-fread(paste0(opt$output,'.',phi_i,'.profiles.chr',i,'.profile'))
			SCORE_temp<-SCORE_temp+profile$SCORESUM
	}
	scores<-cbind(scores, SCORE_temp)

	names(scores)[grepl('SCORE_temp',names(scores))]<-paste0('SCORE_phi',phi_i)

}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

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

###
# Clean up temporary files
###

system(paste0('rm ',opt$output,'.*.profiles.*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
