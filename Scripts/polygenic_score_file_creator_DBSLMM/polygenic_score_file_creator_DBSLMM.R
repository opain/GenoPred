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
make_option("--rscript", action="store", default=NA, type='character',
    help="Path to Rscript binary [optional]"),
make_option("--ld_blocks", action="store", default=NA, type='character',
    help="Path to folder containing LD block information [optional]"),
make_option("--dbslmm", action="store", default=NA, type='character',
    help="Path to DBSLMM directory [optional]"),
make_option("--ldsc", action="store", default=NA, type='character',
    help="Path to LD-score regression binary [required]"),
make_option("--munge_sumstats", action="store", default=NA, type='character',
    help="Path to munge_sumstats.py script [required]"),
make_option("--ldsc_ref", action="store", default=NA, type='character',
    help="Path to LD-score regression reference data 'eur_w_ld_chr' [required]"),
make_option("--hm3_snplist", action="store", default=NA, type='character',
    help="Path to LDSC HapMap3 snplist [required]"),
make_option("--pop_prev", action="store", default=NA, type='numeric',
    help="Population prevelance (if binary) [optional]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--sample_prev", action="store", default=NA, type='numeric', 
    help="Sampling ratio in GWAS [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

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
# polygenic_score_file_creator_DBSLMM.R V1.0
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

system(paste0(opt$munge_sumstats,' --sumstats ',opt$sumstats,' --merge-alleles ',opt$ldsc_ref,'/w_hm3.snplist --ignore BETA,Z --out ', opt$output_dir,'munged_sumstats_temp'))

#####
# Estimate the SNP-heritability
#####

if(opt$pop_prev == 'NA'){
	opt$pop_prev<-NA
}

if(opt$sample_prev == 'NA'){
	opt$sample_prev<-NA
}

if(!is.na(opt$pop_prev) & !is.na(opt$sample_prev)){
  system(paste0(opt$ldsc,' --h2 ',opt$output_dir,'munged_sumstats_temp.sumstats.gz --ref-ld-chr ',opt$ldsc_ref,'/ --w-ld-chr ',opt$ldsc_ref,'/ --out ', opt$output_dir,'ldsc_snp_h2_temp --samp-prev ',opt$sample_prev,' --pop-prev ',opt$pop_prev))

  ldsc_log<-read.table(paste0(opt$output_dir,'ldsc_snp_h2_temp.log'), header=F, sep='&')
  ldsc_h2<-ldsc_log[grepl('Total Liability scale h2', ldsc_log$V1),]
  ldsc_h2<-gsub('Total Liability scale h2: ','', ldsc_h2)
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('SNP-heritability estimate on liability scale = ',ldsc_h2,'.\n',sep='')
  sink()
  
  ldsc_h2<-as.numeric(gsub(' .*','', ldsc_h2))
  
  if(ldsc_h2 > 1){
	ldsc_h2<-1
	
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('SNP-heritability estimate set to 1.\n',sep='')
	sink()
  }
	
} else {
  system(paste0(opt$ldsc,' --h2 ',opt$output_dir,'munged_sumstats_temp.sumstats.gz --ref-ld-chr ',opt$ldsc_ref,'/ --w-ld-chr ',opt$ldsc_ref,'/ --out ', opt$output_dir,'ldsc_snp_h2_temp'))
  
  ldsc_log<-read.table(paste0(opt$output_dir,'ldsc_snp_h2_temp.log'), header=F, sep='&')
  ldsc_h2<-ldsc_log[grepl('Total Observed scale h2', ldsc_log$V1),]
  ldsc_h2<-gsub('Total Observed scale h2: ','', ldsc_h2)
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('SNP-heritability estimate on observed scale = ',ldsc_h2,'.\n',sep='')
  sink()
  
  ldsc_h2<-as.numeric(gsub(' .*','', ldsc_h2))
  
}

if(!is.na(opt$ref_keep)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('ref_keep used to subset reference genotype data.\n')
  sink()
  
  #####
  # Create subset of ref files
  #####
  
  for(i in CHROMS){
      system(paste0(opt$plink,' --bfile ',opt$ref_plink_chr,i,' --keep ',opt$ref_keep,' --make-bed --out ',opt$output_dir,'dbslmm_ref_chr',i))
  }
  
  opt$ref_plink_subset<-paste0(opt$output_dir,'dbslmm_ref_chr')
} else {
  opt$ref_plink_subset<-opt$ref_plink_chr
}

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

# Convert OR to BETA
if(!('BETA' %in% names(GWAS))){
  GWAS$BETA<-log(GWAS$OR)
}

# Rename allele frequency column
if(sum(names(GWAS) == 'FREQ') != 1){
  GWAS$FREQ<-GWAS$REF.FREQ
}

nsnp<-dim(GWAS)[1]
GWAS_N<-mean(GWAS$N)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

# Convert to GEMMA format
if(('CHR' %in% names(GWAS)) & ('BP' %in% names(GWAS))){
  GWAS$N_MISS<-max(GWAS$N)-GWAS$N
  GWAS<-GWAS[,c('CHR','SNP','BP','N_MISS','N','A1','A2','FREQ','BETA','SE','P'),with=F]
  names(GWAS)<-c('chr','rs','ps','n_mis','n_obs','allele1','allele0','af','beta','se','p_wald')
}

# Match allele1 and 0 with A1 and 2 in reference (DBSLMM calls this allele discrepancy)
ref_bim<-NULL
for(i in CHROMS){
  ref_bim<-rbind(ref_bim, fread(paste0(opt$ref_plink_subset, i,'.bim')))
}

GWAS_match<-merge(GWAS, ref_bim[,c('V2','V5','V6'),with=F], by.x=c('rs','allele1','allele0'), by.y=c('V2','V5','V6'))
GWAS_switch<-merge(GWAS, ref_bim[,c('V2','V5','V6'),with=F], by.x=c('rs','allele1','allele0'), by.y=c('V2','V6','V5'))
GWAS_switch$allele_tmp<-GWAS_switch$allele0
GWAS_switch$allele0<-GWAS_switch$allele1
GWAS_switch$allele1<-GWAS_switch$allele_tmp
GWAS_switch$allele_tmp<-NULL
GWAS_switch$beta<--GWAS_switch$beta
GWAS_switch$af<-1-GWAS_switch$af
GWAS<-rbind(GWAS_match, GWAS_switch)

GWAS<-GWAS[order(GWAS$chr, GWAS$ps),]
GWAS<-GWAS[,c('chr','rs','ps','n_mis','n_obs','allele1','allele0','af','beta','se','p_wald'),with=F]

# Write out formatted sumstats
for(i in CHROMS){
  fwrite(GWAS[GWAS$chr == i,], paste0(opt$output_dir,'summary_gemma_chr',i,'.assoc.txt'), sep='\t', col.names=F)
}

rm(GWAS, GWAS_match, GWAS_switch)
gc()

if(!is.na(opt$test)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  test_start.time <- Sys.time()
  cat('Test started at',as.character(test_start.time),'\n')
  sink()
}

#####
# Process sumstats using DBSLMM default mode
#####

system(paste0('chmod 777 ', opt$dbslmm))

for(chr in CHROMS){
  system(paste0(opt$rscript,' ',opt$dbslmm,'/DBSLMM.R --plink ',opt$plink,' --block ',opt$ld_blocks,'/fourier_ls-chr',chr,'.bed --dbslmm ',opt$dbslmm,'/dbslmm --h2 ',ldsc_h2,' --ref ',opt$ref_plink_subset,chr,' --summary ',opt$output_dir,'summary_gemma_chr',chr,'.assoc.txt --n ',round(GWAS_N,0),' --nsnp ',nsnp,' --outPath ',opt$output_dir,' --thread 1'))
}

dbslmm<-list.files(path=opt$output_dir, pattern='.dbslmm.txt')

if(length(dbslmm) != 22){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('At least one chromosome did not complete.\n')
  sink()
}

dbslmm_all<-NULL
for(i in dbslmm){
  dbslmm_all<-rbind(dbslmm_all, fread(paste0(opt$output_dir,'/',i)))
}

write.table(dbslmm_all, paste0(opt$output,'.dbslmm.GW.txt'), quote=F, col.names=F, row.names=F)

if(!is.na(opt$test)){
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  system(paste0('rm ',opt$output_dir,'*.txt'))
  system(paste0('rm ',opt$output_dir,'*.badsnps'))
  system(paste0('rm ',opt$output_dir,'ldsc_snp_h2_temp.log'))
  system(paste0('rm ',opt$output_dir,'munged_sumstats_temp.log'))
  system(paste0('rm ',opt$output_dir,'munged_sumstats_temp.sumstats.gz'))
  if(!is.na(opt$ref_keep)){
    system(paste0('rm ',opt$output_dir,'dbslmm_ref_chr*'))
  }
  q()
}

####
# Calculate mean and sd of polygenic scores
####

# Calculate polygenic scores for reference individuals
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

for(i in CHROMS){
	system(paste0(opt$plink, ' --bfile ',opt$ref_plink_chr,i,' --score ',opt$output,'.dbslmm.GW.txt 1 2 4 sum --out ',opt$output,'.dbslmm.profiles.chr',i,' --memory ',floor(opt$memory*0.7)))
}

# Add up the scores across chromosomes
fam<-fread(paste0(opt$ref_plink_chr,'22.fam'))
scores<-fam[,1:2]
names(scores)<-c('FID','IID')

SCORE_temp<-0
for(i in CHROMS){
		profile<-fread(paste0(opt$output,'.dbslmm.profiles.chr',i,'.profile'))
		SCORE_temp<-SCORE_temp+profile$SCORESUM
}
scores<-cbind(scores, SCORE_temp)

names(scores)[grepl('SCORE_temp',names(scores))]<-'SCORE'

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
if(!is.na(opt$ref_keep)){
  system(paste0('rm ',opt$output_dir,'dbslmm_ref_chr*'))
}
system(paste0('rm ',opt$output_dir,'summary_gemma_chr*.assoc.txt'))
system(paste0('rm ',opt$output_dir,'summary_gemma_chr*.dbslmm.txt'))
system(paste0('rm ',opt$output_dir,'ldsc_snp_h2_temp.log'))
system(paste0('rm ',opt$output_dir,'munged_sumstats_temp.log'))
system(paste0('rm ',opt$output_dir,'munged_sumstats_temp.sumstats.gz'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
