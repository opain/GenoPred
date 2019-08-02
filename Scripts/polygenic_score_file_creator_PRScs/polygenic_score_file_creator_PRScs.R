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
		help="GWAS summary statistics in LDSC format [optional]"),
make_option("--PRScs_path", action="store", default=NA, type='character',
		help="Path to PRScs executable [required]"),
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Number of cores for parallel computing [optional]"),
make_option("--PRScs_ref_path", action="store", default=T, type='character',
		help="Path to PRScs reference [required]"),
make_option("--phi_param", action="store", default='auto', type='character',
		help="Path to PRScs reference [optional]"),
make_option("--python_path", action="store", default=NA, type='character',
		help="Path to python 2.X [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)
system(paste0('mkdir -p ',opt$output_dir))

phi_param<-unlist(strsplit(opt$phi_param,','))

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
GWAS_N<-mean(GWAS$N)
GWAS$N<-NULL
GWAS$P<-2*pnorm(-abs(GWAS$Z))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

GWAS$IUPAC[GWAS$A1 == 'A' & GWAS$A2 =='T' | GWAS$A1 == 'T' & GWAS$A2 =='A']<-'W'
GWAS$IUPAC[GWAS$A1 == 'C' & GWAS$A2 =='G' | GWAS$A1 == 'G' & GWAS$A2 =='C']<-'S'
GWAS$IUPAC[GWAS$A1 == 'A' & GWAS$A2 =='G' | GWAS$A1 == 'G' & GWAS$A2 =='A']<-'R'
GWAS$IUPAC[GWAS$A1 == 'C' & GWAS$A2 =='T' | GWAS$A1 == 'T' & GWAS$A2 =='C']<-'Y'
GWAS$IUPAC[GWAS$A1 == 'G' & GWAS$A2 =='T' | GWAS$A1 == 'T' & GWAS$A2 =='G']<-'K'
GWAS$IUPAC[GWAS$A1 == 'A' & GWAS$A2 =='C' | GWAS$A1 == 'C' & GWAS$A2 =='A']<-'M'

# Extract SNPs that match the reference
GWAS_clean<-NULL
for(i in 1:22){
	bim<-fread(paste0(opt$ref_plink_chr,i,'.bim'))
	
	bim$IUPAC[bim$V5 == 'A' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='A']<-'W'
	bim$IUPAC[bim$V5 == 'C' & bim$V6 =='G' | bim$V5 == 'G' & bim$V6 =='C']<-'S'
	bim$IUPAC[bim$V5 == 'A' & bim$V6 =='G' | bim$V5 == 'G' & bim$V6 =='A']<-'R'
	bim$IUPAC[bim$V5 == 'C' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='C']<-'Y'
	bim$IUPAC[bim$V5 == 'G' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='G']<-'K'
	bim$IUPAC[bim$V5 == 'A' & bim$V6 =='C' | bim$V5 == 'C' & bim$V6 =='A']<-'M'

	bim_GWAS<-merge(bim,GWAS, by.x='V2', by.y='SNP')
	GWAS_clean_temp<-bim_GWAS[bim_GWAS$IUPAC.x == bim_GWAS$IUPAC.y,]
	GWAS_clean<-rbind(GWAS_clean,GWAS_clean_temp)
}

GWAS_clean<-GWAS_clean[,c('V2','A1','A2','Z','P')]
names(GWAS_clean)<-c('SNP','A1','A2','BETA','P')

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After harmonisation with the reference,',dim(GWAS_clean)[1],'variants remain.\n')
sink()

fwrite(GWAS_clean, paste0(opt$output_dir,'GWAS_sumstats_temp.txt'), sep=' ')

rm(GWAS_clean, bim, GWAS_clean_temp, GWAS, bim_GWAS)
gc()

#####
# Process sumstats using PRSsc
#####

# Order chromosomes to be balanced across cores
is.odd <- function(x){x %% 2 != 0}
CHROMS_mat<-matrix(NA,nrow=opt$n_cores, ncol=ceiling(22/opt$n_cores)) 
CHROMS_mat[1:22]<-1:22
for(i in which(is.odd(1:dim(CHROMS_mat)[2]))){CHROMS_mat[,i]<-rev(CHROMS_mat[,i])} 
CHROMS<-as.numeric(CHROMS_mat) 
CHROMS<-CHROMS[!is.na(CHROMS)] 
print(CHROMS)

# Run using PRScs auto, and specifying a range of global shrinkae parameters (1e-6, 1e-4, 1e-2, 1)
foreach(i=CHROMS, .combine=c) %dopar% {
	for(phi_i in phi_param){
		if(phi_i == 'auto'){
			system(paste0('python ',opt$PRScs_path,' --ref_dir=',opt$PRScs_ref_path,' --bim_prefix ',opt$ref_plink_chr,i,' --sst_file=',opt$output_dir,'GWAS_sumstats_temp.txt --n_gwas=',round(GWAS_N,0),' --out_dir=',opt$output,' --chrom=',i))
		} else {
			system(paste0('python ',opt$PRScs_path,' --ref_dir=',opt$PRScs_ref_path,' --bim_prefix ',opt$ref_plink_chr,i,' --phi=',phi_i,' --sst_file=',opt$output_dir,'GWAS_sumstats_temp.txt --n_gwas=',round(GWAS_N,0),' --out_dir=',opt$output,' --chrom=',i))
		}
	}
}

system(paste0('rm ',opt$output_dir,'GWAS_sumstats_temp.txt'))

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
