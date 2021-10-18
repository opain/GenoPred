#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--ref_plink", action="store", default=NA, type='character',
              help="Path to GW reference PLINK files [required]"),
  make_option("--ref_keep", action="store", default=NA, type='character',
              help="Path to keep file for reference [required]"),
  make_option("--ref_pop_scale", action="store", default=NA, type='character',
              help="File containing the population code and location of the keep file [required]"),
  make_option("--plink1", action="store", default='plink', type='character',
              help="Path PLINK v1.9 software binary [required]"),
  make_option("--plink2", action="store", default='plink', type='character',
              help="Path PLINK v2 software binary [required]"),
  make_option("--output", action="store", default='NA', type='character',
              help="Path for output files [required]"),
  make_option("--memory", action="store", default=5000, type='numeric',
              help="Memory limit [optional]"),
  make_option("--sumstats", action="store", default=NA, type='character',
              help="GWAS summary statistics [optional]"),
  make_option("--ldak", action="store", default=NA, type='character',
              help="Path to ldak executable [required]"),
  make_option("--ldak_map", action="store", default=NA, type='character',
              help="Path to ldak map [required]"),
  make_option("--ldak_tag", action="store", default=NA, type='character',
              help="Path to ldak tagging data [required]"),
  make_option("--ldak_highld", action="store", default=NA, type='character',
              help="Path to ldak highld data [required]"),
  make_option("--n_cores", action="store", default=1, type='numeric',
              help="Number of cores for parallel computing [optional]"),
  make_option("--prs_model", action="store", default='mega', type='character',
              help="Model used for deriving SNP-weights [optional]"),
  make_option("--test", action="store", default=NA, type='character',
              help="Specify number of SNPs to include [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

opt$output_dir<-dirname(opt$output)
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
# polygenic_score_file_creator_LDAK.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#####
# Format the sumstats
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in GWAS and harmonising with reference.\n')
sink()

GWAS<-fread(cmd=paste0('zcat ',opt$sumstats), nThread=opt$n_cores)
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

# Calculate Z
if(('BETA' %in% names(GWAS))){
  GWAS$Z<-abs(GWAS$BETA)/GWAS$SE
  GWAS$Z[GWAS$BETA < 0]<- -GWAS$Z[GWAS$BETA < 0]
} else {
  GWAS$Z<-abs(log(GWAS$OR))/GWAS$SE
  GWAS$Z[GWAS$OR < 1]<- -GWAS$Z[GWAS$OR < 1]
}

GWAS$Predictor<-paste0(GWAS$CHR, ':', GWAS$BP)

GWAS<-GWAS[,c('Predictor','A1','A2','N','Z')]
names(GWAS)<-c('Predictor','A1','A2','n','Z')

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

############
# Format reference for LDAK
############

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Formatting reference for LDAK.\n')
sink()

# Must be hg19 with CHR:BP IDs
system(paste0('cp ', opt$ref_plink,'.bed ',opt$output_dir,'/ref.bed'))
system(paste0('cp ', opt$ref_plink,'.fam ',opt$output_dir,'/ref.fam'))
system(paste0("awk < ", opt$ref_plink,".bim '{$2=$1\":\"$4;print $0}' > ", opt$output_dir,'/ref.bim'))

# Insert genetic distances
system(paste0(opt$plink1,' --bfile ',opt$output_dir,'/ref --cm-map ',opt$ldak_map,'/genetic_map_chr@_combined_b37.txt --make-bed --out ',opt$output_dir,'/map'))

system(paste0("cat ", opt$output_dir,"/map.bim | awk '{print $2, $3}' > ", opt$output_dir,"/map.all"))

system(paste0("awk '(NR==FNR){arr[$1]=$2;next}{print $1, $2, arr[$2], $4, $5, $6}' ",opt$output_dir,"/map.all ",opt$output_dir,"/ref.bim > ",opt$output_dir,"/tmp.bim; mv ",opt$output_dir,"/tmp.bim ",opt$output_dir,"/ref.bim"))

system(paste0('rm ',opt$output_dir,'/map*'))

# Extract only SNPs in GWAS
if(!is.na(opt$ref_keep)){
  # Extract ref_keep
  system(paste0(opt$plink1,' --bfile ', opt$output_dir,'/ref --keep ',opt$ref_keep,' --make-bed --out ',opt$output_dir,'/ref_subset --extract ',opt$output_dir,'GWAS_sumstats_temp.txt'))
  system(paste0('mv ',opt$output_dir,'/ref_subset.bed ',opt$output_dir,'/ref.bed'))
  system(paste0('mv ',opt$output_dir,'/ref_subset.bim ',opt$output_dir,'/ref.bim'))
  system(paste0('mv ',opt$output_dir,'/ref_subset.fam ',opt$output_dir,'/ref.fam'))
} else {
  system(paste0(opt$plink1,' --bfile ', opt$output_dir,'/ref --make-bed --out ',opt$output_dir,'/ref_subset --extract ',opt$output_dir,'GWAS_sumstats_temp.txt'))
  system(paste0('mv ',opt$output_dir,'/ref_subset.bed ',opt$output_dir,'/ref.bed'))
  system(paste0('mv ',opt$output_dir,'/ref_subset.bim ',opt$output_dir,'/ref.bim'))
  system(paste0('mv ',opt$output_dir,'/ref_subset.fam ',opt$output_dir,'/ref.fam'))
}

############
# Estimate Per-Predictor Heritabilities
############
# We will use the BLD-LDAK Model, as recommended for human SNP data

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Estimating per-predictor heritabilities.\n')
sink()

# Calculate LDAK weights
system(paste0(opt$ldak,' --cut-weights ',opt$output_dir,'/sections --bfile ',opt$output_dir,'/ref --max-threads ',opt$n_cores))
system(paste0(opt$ldak,' --calc-weights-all ',opt$output_dir,'/sections --bfile ',opt$output_dir,'/ref --max-threads ',opt$n_cores))
system(paste0('mkdir ',opt$output_dir,'/bld'))
system(paste0('cp ',opt$ldak_tag,'/* ',opt$output_dir,'/bld/'))
system(paste0('mv ',opt$output_dir,'/sections/weights.short ',opt$output_dir,'/bld/bld65'))

# Calculate taggings
if(length(CHROMS) != 1){
  system(paste0(opt$ldak,' --calc-tagging ',opt$output_dir,'/bld.ldak --bfile ',opt$output_dir,'/ref --ignore-weights YES --power -.25 --annotation-number 65 --annotation-prefix ',opt$output_dir,'/bld/bld --window-cm 1 --save-matrix YES --max-threads ',opt$n_cores))
} else {
  system(paste0(opt$ldak,' --calc-tagging ',opt$output_dir,'/bld.ldak --bfile ',opt$output_dir,'/ref --ignore-weights YES --power -.25 --annotation-number 65 --annotation-prefix ',opt$output_dir,'/bld/bld --window-cm 1 --chr ',CHROMS,' --save-matrix YES --max-threads ',opt$n_cores))
}

# Calculate Per-Predictor Heritabilities.
system(paste0(opt$ldak,' --sum-hers ',opt$output_dir,'/bld.ldak --tagfile ',opt$output_dir,'/bld.ldak.tagging --summary ',opt$output_dir,'GWAS_sumstats_temp.txt --matrix ',opt$output_dir,'/bld.ldak.matrix --max-threads ',opt$n_cores))

ldak_res_her<-fread(paste0(opt$output_dir,'/bld.ldak.hers'))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('SNP-based heritability estimated to be ',ldak_res_her$Heritability[nrow(ldak_res_her)]," (SD=", ldak_res_her$SD[nrow(ldak_res_her)],").\n",sep='')
sink()

# Identify SNPs in high LD regions
system(paste0(opt$ldak,' --cut-genes ',opt$output_dir,'/highld --bfile ',opt$output_dir,'/ref --genefile ',opt$ldak_highld,' --max-threads ',opt$n_cores))

###################
# Run using full reference.
###################

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Running using full reference.\n')
sink()

######
# Calculate predictor-predictor correlations
######

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating predictor-predictor correlations.\n')
sink()

if(length(CHROMS) !=1){
  for(chr in CHROMS){
    system(paste0(opt$ldak,' --calc-cors ',opt$output_dir,'/cors_full',chr,' --bfile ',opt$output_dir,'/ref --window-cm 3 --chr ',chr,' --max-threads ',opt$n_cores))
  }
  
  write.table(paste0(opt$output_dir,'/cors_full',CHROMS), paste0(opt$output_dir,'/cors_full_list.txt'), col.names=F, row.names=F, quote=F)
  
  system(paste0(opt$ldak,' --join-cors ',opt$output_dir,'/cors_full --corslist ',opt$output_dir,'/cors_full_list.txt --max-threads ',opt$n_cores))
} else {
  system(paste0(opt$ldak,' --calc-cors ',opt$output_dir,'/cors_full --bfile ',opt$output_dir,'/ref --window-cm 3 --chr ',CHROMS,' --max-threads ',opt$n_cores))
}

######
# Estimate effect sizes for training and full prediction models.
######

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Running MegaPRS: ',opt$prs_model,' model.\n', sep='')
sink()

system(paste0(opt$ldak,' --mega-prs ',opt$output_dir,'/mega_full --model ',opt$prs_model,' --bfile ',opt$output_dir,'/ref --cors ',opt$output_dir,'/cors_full --ind-hers ',opt$output_dir,'/bld.ldak.ind.hers --summary ',opt$output_dir,'GWAS_sumstats_temp.txt --one-sums YES --window-cm 1 --allow-ambiguous YES --max-threads ',opt$n_cores))

# Save the parameters file
system(paste0('cp ',opt$output_dir,'/mega_full.parameters ',opt$output,'.model_param.txt'))

# Sum of per SNP heritability is different from SNP-heritability, due to removal of variants with non-positive heritability

################
# Run using subset reference for pseudovalidation
################

############
# Create pseudo summaries
############

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Creating pseudosummaries.\n')
sink()

# Split reference into three
system(paste0("awk < ",opt$output_dir,"/ref.fam '(NR%3==1){print $0 > \"",opt$output_dir,"/keepa\"}(NR%3==2){print $0 > \"",opt$output_dir,"/keepb\"}(NR%3==0){print $0 > \"",opt$output_dir,"/keepc\"}'"))

# Create pseudo summaries
system(paste0(opt$ldak,' --pseudo-summaries ',opt$output_dir,'GWAS_sumstats_temp.pseudo --bfile ',opt$output_dir,'/ref --summary ',opt$output_dir,'GWAS_sumstats_temp.txt --training-proportion .9 --keep ',opt$output_dir,'/keepa --allow-ambiguous YES --max-threads ',opt$n_cores))

######
# Calculate predictor-predictor correlations (can be split by chromosome).
######

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating predictor-predictor correlations.\n')
sink()

if(length(CHROMS) !=1){
  for(chr in CHROMS){
    system(paste0(opt$ldak,' --calc-cors ',opt$output_dir,'/cors_subset',chr,' --bfile ',opt$output_dir,'/ref --window-cm 3 --keep ',opt$output_dir,'/keepb --chr ',chr,' --max-threads ',opt$n_cores))
  }
  
  write.table(paste0(opt$output_dir,'/cors_subset',CHROMS), paste0(opt$output_dir,'/cors_subset_list.txt'), col.names=F, row.names=F, quote=F)
  
  system(paste0(opt$ldak,' --join-cors ',opt$output_dir,'/cors_subset --corslist ',opt$output_dir,'/cors_subset_list.txt --max-threads ',opt$n_cores))
} else {
  system(paste0(opt$ldak,' --calc-cors ',opt$output_dir,'/cors_subset --bfile ',opt$output_dir,'/ref --window-cm 3 --keep ',opt$output_dir,'/keepb --chr ',CHROMS,' --max-threads ',opt$n_cores))
}

######
# Estimate effect sizes for training and full prediction models.
######

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Running MegaPRS: ',opt$prs_model,' model.\n', sep='')
sink()

system(paste0(opt$ldak,' --mega-prs ',opt$output_dir,'/mega_subset --model ',opt$prs_model,' --bfile ',opt$output_dir,'/ref --cors ',opt$output_dir,'/cors_subset --ind-hers ',opt$output_dir,'/bld.ldak.ind.hers --summary ',opt$output_dir,'GWAS_sumstats_temp.pseudo.train.summaries --one-sums YES --window-cm 1 --allow-ambiguous YES --max-threads ',opt$n_cores))

# Sum of per SNP heritability is different from SNP-heritability, due to removal of variants with non-positive heritability

######
# Test training prediction models
######
# Equivalent to pseudovalidation?

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Running pseudovalidation.\n')
sink()

if(file.exists(paste0(opt$output_dir,'/highld/genes.predictors.used'))){
  system(paste0(opt$ldak,' --calc-scores ',opt$output_dir,'/mega_subset --bfile ',opt$output_dir,'/ref --scorefile ',opt$output_dir,'/mega_subset.effects --summary ',opt$output_dir,'GWAS_sumstats_temp.pseudo.test.summaries --power 0 --final-effects ',opt$output_dir,'/mega_subset.effects --keep ',opt$output_dir,'/keepc --allow-ambiguous YES --exclude ',opt$output_dir,'/highld/genes.predictors.used --max-threads ',opt$n_cores))
} else {
  system(paste0(opt$ldak,' --calc-scores ',opt$output_dir,'/mega_subset --bfile ',opt$output_dir,'/ref --scorefile ',opt$output_dir,'/mega_subset.effects --summary ',opt$output_dir,'GWAS_sumstats_temp.pseudo.test.summaries --power 0 --final-effects ',opt$output_dir,'/mega_subset.effects --keep ',opt$output_dir,'/keepc --allow-ambiguous YES --max-threads ',opt$n_cores))
}

# Identify the best fitting model
ldak_res_cors<-fread(paste0(opt$output_dir,'/mega_subset.cors'), nThread=opt$n_cores)
best_score<-ldak_res_cors[ldak_res_cors$V2 == max(ldak_res_cors$V2),]

# Save the pseudovalidation results
system(paste0('cp ',opt$output_dir,'/mega_subset.cors ',opt$output,'.pseudoval.txt'))
sink(file = paste(opt$output,'.log',sep=''), append = T)

cat('Model ',gsub('Score_','',best_score$V1[1]),' is identified as the best with correlation of ',best_score$V2,'.\n', sep='')
sink()

if(!is.na(opt$test)){
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  system(paste0('rm ',opt$output_dir,'cors*'))
  system(paste0('rm -r ',opt$output_dir,'bld'))
  system(paste0('rm ',opt$output_dir,'bld*'))
  system(paste0('rm ',opt$output_dir,'GWAS_sums*'))
  system(paste0('rm ',opt$output_dir,'keepa'))
  system(paste0('rm ',opt$output_dir,'keepb'))
  system(paste0('rm ',opt$output_dir,'keepc'))
  system(paste0('rm -r ',opt$output_dir,'highld'))
  system(paste0('rm -r ',opt$output_dir,'sections'))
  system(paste0('rm ',opt$output_dir,'ref*'))
  system(paste0('rm ',opt$output_dir,'mega*'))
  q()
}

######
# Format final score file
######

# Read in the scores
score<-fread(paste0(opt$output_dir,'/mega_full.effects'), nThread=opt$n_cores)

# Change IDs to RSIDs
bim<-fread(paste0(opt$ref_plink,'.bim'), nThread=opt$n_cores)
bim$Predictor<-paste0(bim$V1,':',bim$V4)
score<-merge(score, bim[,c('Predictor','V2'),with=F], by='Predictor')
score<-score[,c('V2','A1',names(score)[grepl('Model', names(score))]), with=F]
names(score)[1]<-'SNP'
names(score)[grepl('Model', names(score))]<-paste0('SCORE_ldak_',names(score)[grepl('Model', names(score))])

fwrite(score, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

if(file.exists(paste0(opt$output,'.score.gz'))){
  system(paste0('rm ',opt$output,'.score.gz'))
}

system(paste0('gzip ',opt$output,'.score'))

####
# Calculate mean and sd of polygenic scores
####

# Calculate polygenic scores for reference individuals
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

system(paste0(opt$plink2, ' --bfile ',opt$ref_plink,' --score ',opt$output,'.score.gz header-read --score-col-nums 3-',ncol(score),' --out ',opt$output,'.profiles --threads ',opt$n_cores,' --memory ',floor(opt$memory*0.7)))

# Add up the scores across chromosomes
sscore<-fread(paste0(opt$output,'.profiles.sscore'), nThread=opt$n_cores)
scores<-sscore[,grepl('SCORE_', names(sscore)),with=F]
scores<-as.matrix(scores*sscore$NMISS_ALLELE_CT)

scores<-data.table(sscore[,1:2,with=F],
                   scores)

names(scores)<-c('FID','IID',names(score)[-1:-2])

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

# Calculate the mean and sd of scores for each population specified in pop_scale
pop_keep_files<-read.table(opt$ref_pop_scale, header=F, stringsAsFactors=F)

for(k in 1:dim(pop_keep_files)[1]){
  pop<-pop_keep_files$V1[k]
  keep<-fread(pop_keep_files$V2[k], header=F, nThread=opt$n_cores)
  scores_keep<-scores[(scores$FID %in% keep$V1),]
  
  ref_scale<-data.frame(	Param=names(scores_keep[,-1:-2]),
                         Mean=round(sapply(scores_keep[,-1:-2], function(x) mean(x)),3),
                         SD=round(sapply(scores_keep[,-1:-2], function(x) sd(x)),3))
  
  fwrite(ref_scale, paste0(opt$output,'.',pop,'.scale'), sep=' ')
}

###
# Clean up temporary files
###

system(paste0('rm ',opt$output,'*.profiles.*'))
system(paste0('rm ',opt$output_dir,'cors*'))
system(paste0('rm -r ',opt$output_dir,'bld'))
system(paste0('rm ',opt$output_dir,'bld*'))
system(paste0('rm ',opt$output_dir,'GWAS_sums*'))
system(paste0('rm ',opt$output_dir,'keepa'))
system(paste0('rm ',opt$output_dir,'keepb'))
system(paste0('rm ',opt$output_dir,'keepc'))
system(paste0('rm -r ',opt$output_dir,'highld'))
system(paste0('rm -r ',opt$output_dir,'sections'))
system(paste0('rm ',opt$output_dir,'ref*'))
system(paste0('rm ',opt$output_dir,'mega*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()

