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
make_option("--sumstats", action="store", default=NA, type='character',
		help="GWAS summary statistics in LDSC format [required]"),
make_option("--gcta", action="store", default=NA, type='character',
		help="Path to GCTA binary [required]"),
make_option("--ldsc", action="store", default=NA, type='character',
    help="Path to LD-score regression binary [required]"),
make_option("--munge_sumstats", action="store", default=NA, type='character',
    help="Path to munge_sumstats.py script [required]"),
make_option("--ldsc_ref", action="store", default=NA, type='character',
    help="Path to LD-score regression reference data 'eur_w_ld_chr' [required]"),
make_option("--hm3_snplist", action="store", default=NA, type='character',
    help="Path to LDSC HapMap3 snplist [required]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--prune_hla", action="store", default=T, type='logical',
		help="Retain only top assocaited variant in HLA region [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)
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
# polygenic_score_file_creator_SBLUP.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#####
# Munge_sumstats
#####

system(paste0(opt$munge_sumstats,' --sumstats ',opt$sumstats,' --merge-alleles ',opt$ldsc_ref,'/w_hm3.snplist --out ', opt$output_dir,'munged_sumstats_temp'))

#####
# Estimate the SNP-heritability
#####

system(paste0(opt$ldsc,' --h2 ',opt$output_dir,'munged_sumstats_temp.sumstats.gz --ref-ld-chr ',opt$ldsc_ref,'/ --w-ld-chr ',opt$ldsc_ref,'/ --out ', opt$output_dir,'ldsc_snp_h2_temp'))

ldsc_log<-read.table(paste0(opt$output_dir,'ldsc_snp_h2_temp.log'), header=F, sep='&')
ldsc_h2<-ldsc_log[grepl('Total Observed scale h2', ldsc_log$V1),]
ldsc_h2<-gsub('Total Observed scale h2: ','', ldsc_h2)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('SNP-heritability estimate = ',ldsc_h2,'.\n',sep='')
sink()

ldsc_h2<-as.numeric(gsub(' .*','', ldsc_h2))

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

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

nsnp<-dim(GWAS)[1]

###
# Change to COJO format
###

# If OR present, calculate BETA
if(sum(names(GWAS) == 'OR') == 1){
  GWAS$BETA<-log(GWAS$OR)
}

# Rename allele frequency column
if(sum(names(GWAS) == 'FREQ') == 1){
  GWAS$MAF<-GWAS$FREQ
} else {
  GWAS$MAF<-GWAS$REF.FREQ
}

GWAS<-GWAS[,c('SNP','A1','A2','MAF','BETA','SE','P','N'),with=F]
names(GWAS)<-c('SNP','A1','A2','freq','b','se','p','N')

fwrite(GWAS, paste0(opt$output_dir,'GWAS_sumstats_COJO.txt'), sep=' ', na = "NA", quote=F)

if(!is.na(opt$test)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  test_start.time <- Sys.time()
  cat('Test started at',as.character(test_start.time),'\n')
  sink()
}

#####
# Run GCTA SBLUP
#####

system(paste0(opt$gcta,' --bfile ',opt$ref_plink,' --maf 0.01 --keep ',opt$ref_keep,' --cojo-file ',opt$output_dir,'GWAS_sumstats_COJO.txt --cojo-sblup ',nsnp*(1/ldsc_h2-1),' --cojo-wind 1000 --thread-num ',opt$n_cores,' --out ',opt$output_dir,'GWAS_sumstats_SBLUP'))

if(!is.na(opt$test)){
  if(!file.exists(paste0(opt$output_dir,'GWAS_sumstats_SBLUP.sblup.cojo'))){
    q()
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  system(paste0('rm ',opt$output_dir,'ldsc_snp_h2_temp.log'))
  system(paste0('rm ',opt$output_dir,'munged_sumstats_temp.log'))
  system(paste0('rm ',opt$output_dir,'munged_sumstats_temp.sumstats.gz'))
  system(paste0('rm ',opt$output_dir,'GWAS_sumstats_COJO.txt'))
  system(paste0('rm ',opt$output_dir,'GWAS_sumstats_SBLUP.log'))
  system(paste0('rm ',opt$output_dir,'GWAS_sumstats_SBLUP.sblup.cojo'))
  q()
}

####
# Calculate mean and sd of polygenic scores at each threshold
####

# Calculate polygenic scores for reference individuals
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

system(paste0(opt$plink, ' --bfile ',opt$ref_plink,' --score ',opt$output_dir,'GWAS_sumstats_SBLUP.sblup.cojo 1 2 4 sum --out ',opt$output_dir,'ref.profiles --memory ',floor(opt$memory*0.7)))

# Read in the reference scores
scores<-fread(paste0(opt$output_dir,'ref.profiles.profile'))

# Calculate the mean and sd of scores for each population specified in pop_scale
pop_keep_files<-read.table(opt$ref_pop_scale, header=F, stringsAsFactors=F)

for(k in 1:dim(pop_keep_files)[1]){
	pop<-pop_keep_files$V1[k]
	keep<-fread(pop_keep_files$V2[k], header=F)
	scores_keep<-scores[(scores$FID %in% keep$V1),]

	ref_scale<-data.frame(	Mean=round(mean(scores_keep$SCORESUM),3),
							SD=round(sd(scores_keep$SCORESUM),3))

	fwrite(ref_scale, paste0(opt$output,'.',pop,'.scale'), sep=' ')
}

###
# Clean up temporary files
###

system(paste0('rm ',opt$output_dir,'ref.profiles.*'))
system(paste0('rm ',opt$output_dir,'ldsc_snp_h2_temp.log'))
system(paste0('rm ',opt$output_dir,'munged_sumstats_temp.log'))
system(paste0('rm ',opt$output_dir,'munged_sumstats_temp.sumstats.gz'))
system(paste0('rm ',opt$output_dir,'GWAS_sumstats_COJO.txt'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
