#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--ref_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK files [required]"),
make_option("--ref_keep", action="store", default=NA, type='character',
		help="Keep file to subset individuals in reference for clumping [required]"),
make_option("--ref_freq_chr", action="store", default=NA, type='character',
    help="Path to per chromosome PLINK frequency (.frq) files [optional]"),
make_option("--ref_maf", action="store", default=NA, type='numeric',
    help="Minor allele frequency threshold to be applied based on ref_freq_chr [optional]"),
make_option("--ref_pop_scale", action="store", default=NA, type='character',
		help="File containing the population code and location of the keep file [required]"),
make_option("--plink", action="store", default='plink', type='character',
		help="Path PLINK software binary [required]"),
make_option("--prsice_path", action="store", default=NA, type='character',
    help="Path to PRSice. [optional]"),
make_option("--rscript", action="store", default='Rscript', type='character',
    help="Path to Rscript [optional]"),
make_option("--output", action="store", default='./Output', type='character',
		help="Path for output files [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]"),
make_option("--sumstats", action="store", default=NA, type='character',
		help="GWAS summary statistics in LDSC format [required]"),
make_option("--covar", action="store", default=NA, type='character',
		help="File containing covariates to be regressed from the scores [optional]"),
make_option("--pTs", action="store", default='1e-8,1e-6,1e-4,1e-2,0.1,0.2,0.3,0.4,0.5,1', type='character',
		help="List of p-value thresholds for scoring [optional]"),
make_option("--dense", action="store", default=F, type='logical',
    help="Specify as T for dense thresholding. pTs then interpretted as seq() command wih default 5e-8,1,5e-4 [optional]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--prune_hla", action="store", default=T, type='logical',
		help="Retain only top assocaited variant in HLA region [optional]")
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
# polygenic_score_file_creator_dense.R V1.0
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
cat('Reading in GWAS.\n')
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
    for(i in CHROMS){
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

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

GWAS<-GWAS[,c('SNP','A1','A2','BETA','P')]

fwrite(GWAS, paste0(opt$output_dir,'GWAS_sumstats_temp.txt'), sep=' ')

if(opt$prune_hla == T){
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Extracted top variant in HLA region.\n')
	sink()

	bim<-fread(paste0(opt$ref_plink_chr,'6.bim'))
	bim_GWAS<-merge(bim,GWAS, by.x='V2',by.y='SNP')
	bim_GWAS_hla<-bim_GWAS[bim_GWAS$V4 > 28e6 & bim_GWAS$V4 < 34e6,]
	bim_GWAS_hla_excl<-bim_GWAS_hla[bim_GWAS_hla$P != min(bim_GWAS_hla$P),]
	write.table(bim_GWAS_hla_excl, paste0(opt$output_dir,'hla_exclude.txt'), col.names=F, row.names=F, quote=F)
}

if(!is.na(opt$test)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  test_start.time <- Sys.time()
  cat('Test started at',as.character(test_start.time),'\n')
  sink()
}

#####
# Clump SNPs in GWAS based on LD in the reference
#####

if(!is.na(opt$ref_keep)){ 
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Keeping individuals in',opt$ref_keep,'for clumping.\n')
	sink()
} else {
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('Using full reference for clumping.\n')
	sink()
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Clumping GWAS based on the reference...')
sink()

if(!is.na(opt$ref_keep)){ 
	if(opt$prune_hla == F){
		for(i in CHROMS){
			system(paste0(opt$plink,' --bfile ',opt$ref_plink_chr,i,' --keep ',opt$ref_keep,' --clump ',opt$output_dir,'GWAS_sumstats_temp.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 250 --out ',opt$output_dir,'GWAS_sumstats_temp_clumped_chr',i,'.txt --memory ',floor(opt$memory*0.7)))
		}
	} else {
		for(i in CHROMS){
			system(paste0(opt$plink,' --bfile ',opt$ref_plink_chr,i,' --keep ',opt$ref_keep,' --exclude ',opt$output_dir,'hla_exclude.txt --clump ',opt$output_dir,'GWAS_sumstats_temp.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 250 --out ',opt$output_dir,'GWAS_sumstats_temp_clumped_chr',i,'.txt --memory ',floor(opt$memory*0.7)))
		}
	}
} else {
	if(opt$prune_hla == F){
		for(i in CHROMS){
			system(paste0(opt$plink,' --bfile ',opt$ref_plink_chr,i,' --clump ',opt$output_dir,'GWAS_sumstats_temp.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 250 --out ',opt$output_dir,'GWAS_sumstats_temp_clumped_chr',i,'.txt --memory ',floor(opt$memory*0.7)))
		}
	} else {
		for(i in CHROMS){
			system(paste0(opt$plink,' --bfile ',opt$ref_plink_chr,i,' --exclude ',opt$output_dir,'hla_exclude.txt --clump ',opt$output_dir,'GWAS_sumstats_temp.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 250 --out ',opt$output_dir,'GWAS_sumstats_temp_clumped_chr',i,'.txt --memory ',floor(opt$memory*0.7)))
		}
	}
}

GWAS_clumped_all<-NULL
for(i in CHROMS){
  clumped<-fread(paste0(opt$output_dir,'GWAS_sumstats_temp_clumped_chr',i,'.txt.clumped'))
  clumped_SNPs<-clumped$SNP
  GWAS_clumped_temp<-GWAS[(GWAS$SNP %in% clumped_SNPs),]
  GWAS_clumped_all<-rbind(GWAS_clumped_all, GWAS_clumped_temp)
}

write.table(GWAS_clumped_all, paste0(opt$output,'.GWAS_sumstats_clumped.txt'), col.names=T, row.names=F, quote=F)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Clean up temporary files
###

system(paste0('rm ',opt$output_dir,'GWAS*'))
if(opt$prune_hla == T){
	system(paste0('rm ',opt$output_dir,'hla_exclude.txt'))
}

####
# Calculate mean and sd of polygenic scores at each threshold
####

# Calculate polygenic scores for reference individuals
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

system(paste0(opt$rscript,' ',opt$prsice,'/PRSice.R --prsice ',opt$prsice,'/PRSice_linux --base ',opt$output,'.GWAS_sumstats_clumped.txt --target ',opt$ref_plink_chr,"# --thread 1 --lower 1e-8 --stat BETA --score sum --binary-target F --no-clump --no-regress --out ",opt$output,'ref_score'))

if(!is.na(opt$test)){
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  
  system(paste0('rm ',opt$output,'ref_score.all.score'))
  system(paste0('rm ',opt$output,'ref_score.log'))
  system(paste0('rm ',opt$output,'ref_score.prsice'))
  system(paste0('rm ',opt$output_dir,'hla_exclude.txt'))
  system(paste0('rm ',opt$output_dir,'GWAS_sumstats_temp.txt'))
  system(paste0('rm ',opt$output,'.GWAS_sumstats_clumped.txt'))
  q()
}

# Read in the scores
if(file.exists(paste0(opt$output,'ref_score.all.score'))){
  scores<-fread(paste0(opt$output,'ref_score.all.score'), header=T)
} else {
  scores<-fread(paste0(opt$output,'ref_score.all_score'), header=T)
}

names(scores)[-1:-2]<-paste0('SCORE_',names(scores)[-1:-2])

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

# Calculate the mean and sd of scores for each population specified in pop_scale
pop_keep_files<-read.table(opt$ref_pop_scale, header=F, stringsAsFactors=F)

for(k in 1:dim(pop_keep_files)[1]){
	pop<-pop_keep_files$V1[k]
	keep<-fread(pop_keep_files$V2[k], header=F)
	scores_keep<-scores[(scores$FID %in% keep$V1),]

	ref_scale<-data.frame(	pT=names(scores_keep[,-1:-2]),
													Mean=round(sapply(scores_keep[,-1:-2], function(x) mean(x)),3),
													SD=round(sapply(scores_keep[,-1:-2], function(x) sd(x)),3))

	fwrite(ref_scale, paste0(opt$output,'.',pop,'.scale'), sep=' ')
}

###
# Clean up temporary files
###

if(file.exists(paste0(opt$output,'ref_score.all.score'))){
  system(paste0('rm ',opt$output,'ref_score.all.score'))
} else {
  system(paste0('rm ',opt$output,'ref_score.all_score'))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
