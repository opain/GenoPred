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
make_option("--output", action="store", default='./Output', type='character',
		help="Path for output files [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]"),
make_option("--sumstats", action="store", default=NA, type='character',
		help="GWAS summary statistics in LDSC format [required]"),
make_option("--covar", action="store", default=NA, type='character',
    help="File containing covariates to be regressed from the scores [optional]"),
make_option("--annot", action="store", default=NA, type='character',
    help="File containing SNP annotations [optional]"),
make_option("--annot_col_name", action="store", default=NA, type='character',
    help="Column name for annotations [optional]"),
make_option("--annot_weight", action="store", default=NA, type='numeric',
    help="Weight to be assigned to annotation [optional]"),
make_option("--pTs", action="store", default='1e-8,1e-6,1e-4,1e-2,0.1,0.2,0.3,0.4,0.5,1', type='character',
		help="List of p-value thresholds for scoring [optional]"),
make_option("--extract", action="store", default=NA, type='character',
    help="File listing SNPs to extract for polygenic scoring [optional]"),
make_option("--nested", action="store", default=T, type='logical',
    help="Specify as F to use non-overlapping p-value intervals [optional]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--prune_hla", action="store", default=T, type='logical',
		help="Retain only top assocaited variant in HLA region [optional]")
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
# polygenic_score_file_creator_weighted.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#####
# Format pT option
#####

opt$pTs<-as.numeric(unlist(strsplit(opt$pTs,',')))

#####
# Read in sumstats
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

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

if(!is.na(opt$extract)){
	extract_snplist<-fread(opt$extract, header=F)$V1
	GWAS<-GWAS[(GWAS$SNP %in% extract_snplist),]
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('After applying the extract file,',dim(GWAS)[1],'variants remain.\n')
	sink()
}

###
# Adjust BETA and P for annotation
###

annot<-fread(opt$annot)
if(!is.na(opt$annot_col_name)){
  names(annot)[names(annot) == opt$annot_col_name]<-'annot'
  annot<-annot[,c('SNP','annot'),with=F]
}
GWAS<-merge(GWAS, annot, by='SNP')
GWAS$Z<-abs(qnorm(GWAS$P/2))
GWAS$SE<-GWAS$BETA/GWAS$Z
GWAS$BETA<-GWAS$BETA*(GWAS$annot^opt$annot_weight)
GWAS$Z<-GWAS$BETA/GWAS$SE
GWAS$P<-2*pnorm(-abs(GWAS$Z))

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

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

#####
# Create score files
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Creating score files...')
sink()

GWAS_clumped_all<-NULL
for(i in CHROMS){
	clumped<-fread(paste0(opt$output_dir,'GWAS_sumstats_temp_clumped_chr',i,'.txt.clumped'))
	clumped_SNPs<-clumped$SNP
	GWAS_clumped_temp<-GWAS[(GWAS$SNP %in% clumped_SNPs),]
	fwrite(GWAS_clumped_temp[,c('SNP','A1','BETA')], paste0(opt$output,'.chr',i,'.score'), sep=' ', col.names=F)
	fwrite(GWAS_clumped_temp[,c('SNP','P')], paste0(opt$output,'.chr',i,'.range_values'), sep=' ', col.names=F)
	GWAS_clumped_all<-rbind(GWAS_clumped_all, GWAS_clumped_temp[,c('SNP','P')])
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

#####
# Create range list file
#####
if(opt$nested == T){
  range_list<-data.frame(	Name=paste0('S',1:length(opt$pTs)),
                          pT0=0,
                          pT1=opt$pTs)
} else {
  range_list<-data.frame(	Name=paste0('S',1:length(opt$pTs)),
                          pT0=c(0,opt$pTs[-length(opt$pTs)]),
                          pT1=opt$pTs)
}

if(!is.na(opt$test)){
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  system(paste0('rm ',opt$output,'.chr*'))
  system(paste0('rm ',opt$output_dir,'GWAS_sumstats_temp_clumped_chr*'))
  system(paste0('rm ',opt$output_dir,'hla_exclude.txt'))
  system(paste0('rm ',opt$output_dir,'GWAS_sumstats_temp.txt'))
  q()
}

###
# Count the number of SNPs to be included at each pT
###

for(i in 1:dim(range_list)[1]){
	range_list$NSNP[i]<-sum(GWAS_clumped_all$P > range_list$pT0[i] & GWAS_clumped_all$P < range_list$pT1[i])
}

fwrite(range_list, paste0(opt$output,'.NSNP_per_pT'),sep='\t')

# Output range list file with ranges containing at least 1 SNP
range_list<-range_list[range_list$NSNP > 0,]	
range_list$NSNP<-NULL
fwrite(range_list, paste0(opt$output,'.range_list'), col.names=F, sep=' ')

rm(GWAS_clumped_all)
gc()

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

for(i in CHROMS){
	system(paste0(opt$plink, ' --bfile ',opt$ref_plink_chr,i,' --score ',opt$output,'.chr',i,'.score sum --q-score-range ',opt$output,'.range_list ',opt$output,'.chr',i,'.range_values --out ',opt$output,'.profiles.chr',i,' --memory ',floor(opt$memory*0.7)))
}

# Add up the scores across chromosomes
fam<-fread(paste0(opt$ref_plink_chr,'22.fam'))
scores<-fam[,1:2]
names(scores)<-c('FID','IID')

for(k in 1:dim(range_list)[1]){
SCORE_temp<-0
	for(i in CHROMS){
		if(file.exists(paste0(opt$output,'.profiles.chr',i,'.',range_list$Name[k],'.profile'))){			
			profile<-fread(paste0(opt$output,'.profiles.chr',i,'.',range_list$Name[k],'.profile'))
			SCORE_temp<-SCORE_temp+profile$SCORE
		}
	}
	scores<-cbind(scores, SCORE_temp)
	names(scores)[k+2]<-paste0('SCORE_',range_list$pT1[k])
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

# Regress out covariates if specified
if(!is.na(opt$covar) == T){
	covar<-fread(opt$covar)
	scores_covar<-merge(scores,covar,by=c('FID','IID'))
	
	# Scale the covariates so coeficients correspond to covariates in target samples
	for(i in names(scores_covar)[grepl('PC',names(scores_covar))]){
		scores_covar[[i]]<-as.numeric(scale(scores_covar[[i]]))
	}
	
	models<-list()
	scores_resid<-data.frame(scores_covar[,c('FID','IID')])
	for(i in names(scores[,-1:-2])){
		models[[i]]<-lm(as.formula(paste0('scores_covar[[i]] ~ ',paste(names(scores_covar)[grepl('PC', names(scores_covar))], collapse=' + '))), data=scores_covar)
		scores_resid[[i]]<-resid(models[[i]])
	}
	
saveRDS(models, paste0(opt$output,'.models.rds'))
scores<-scores_resid

}

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

system(paste0('rm ',opt$output,'.profiles.*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
