#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
	make_option("--ref_plink_gw", action="store", default=NA, type='character',
			help="Path to genome-wide reference PLINK files [required]"),
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
			help="GWAS summary statistics in LDSC format [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(lassosum)
setwd(system.file("data", package="lassosum"))

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# polygenic_score_file_creator_lassosum.R V1.0
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
GWAS$P<-2*pnorm(-abs(GWAS$Z))
GWAS$N<-NULL

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

GWAS$IUPAC[GWAS$A1 == 'A' & GWAS$A2 =='T' | GWAS$A1 == 'T' & GWAS$A2 =='A']<-'W'
GWAS$IUPAC[GWAS$A1 == 'C' & GWAS$A2 =='G' | GWAS$A1 == 'G' & GWAS$A2 =='C']<-'S'
GWAS$IUPAC[GWAS$A1 == 'A' & GWAS$A2 =='G' | GWAS$A1 == 'G' & GWAS$A2 =='A']<-'R'
GWAS$IUPAC[GWAS$A1 == 'C' & GWAS$A2 =='T' | GWAS$A1 == 'T' & GWAS$A2 =='C']<-'Y'
GWAS$IUPAC[GWAS$A1 == 'G' & GWAS$A2 =='T' | GWAS$A1 == 'T' & GWAS$A2 =='G']<-'K'
GWAS$IUPAC[GWAS$A1 == 'A' & GWAS$A2 =='C' | GWAS$A1 == 'C' & GWAS$A2 =='A']<-'M'

#####
# Extract SNPs that match the reference
#####
bim<-fread(paste0(opt$ref_plink_gw,'.bim'), nThread=1)

bim$IUPAC[bim$V5 == 'A' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='A']<-'W'
bim$IUPAC[bim$V5 == 'C' & bim$V6 =='G' | bim$V5 == 'G' & bim$V6 =='C']<-'S'
bim$IUPAC[bim$V5 == 'A' & bim$V6 =='G' | bim$V5 == 'G' & bim$V6 =='A']<-'R'
bim$IUPAC[bim$V5 == 'C' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='C']<-'Y'
bim$IUPAC[bim$V5 == 'G' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='G']<-'K'
bim$IUPAC[bim$V5 == 'A' & bim$V6 =='C' | bim$V5 == 'C' & bim$V6 =='A']<-'M'

bim_GWAS<-merge(bim, GWAS, by.x='V2', by.y='SNP')
GWAS_clean<-bim_GWAS[bim_GWAS$IUPAC.x == bim_GWAS$IUPAC.y,]

GWAS_clean<-GWAS_clean[,c('V1','V2','V4','A1','A2','Z','P')]
names(GWAS_clean)<-c('CHR','SNP','BP','A1','A2','Z','P')

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After harmonisation with the reference,',dim(GWAS_clean)[1],'variants remain.\n')
sink()

rm(bim_GWAS, bim, GWAS)
gc()

#####
# Calculate correlation between SNP and phenotype 
#####
cor <- p2cor(p = GWAS_clean$P, n = GWAS_N, sign=GWAS_clean$Z)

#####
# Perform lassosum to shrink effects using a range of parameters
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Running lassosum pipeline...')
sink()

out<-lassosum.pipeline(cor=cor, chr=GWAS_clean$CHR, pos=GWAS_clean$BP, 
                       A1=GWAS_clean$A1, A2=GWAS_clean$A2,
                       ref.bfile=opt$ref_plink_gw, keep.ref=opt$ref_keep, 
                       LDblocks = 'EUR.hg19')

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

rm(cor)
gc()

#####
# Perform pseudovalidation to idenitfy the best p-value threshold
#####
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Idenitfying best parameters via pseudovalidation...')
sink()

png(paste0(opt$output,'.pseudovalidate.png'), unit='px', res=300, height=2000, width=2000)
v <- pseudovalidate(out)
dev.off()

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

# Subset the validated lassosum model
out2 <- subset(out, s=v$best.s, lambda=v$best.lambda)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Pseudovalidated parameters:
s = ',out2$s,'
lambda = ',out2$lambda,'\n',sep='')
sink()

rm(v)
gc()

# Write out a score file
score_file<-data.frame(SNP=GWAS_clean$SNP[out$sumstats$order],
                       out$sumstats[c('chr','pos','A1','A2')])
score_file<-score_file[,c('chr','SNP','pos','A1','A2')]
names(score_file)<-c('CHR','SNP','pos','A1','A2')

for(i in 1:length(out$s)){
  for(k in 1:length(out$lambda)){
    score_file_tmp<-data.frame(out$beta[[i]][,k])
    names(score_file_tmp)<-paste0('s',out$s[i],'_lambda',out$lambda[k])
    score_file_tmp<-cbind(score_file,score_file_tmp)
    write.table(score_file_tmp, paste0(opt$output,'.s',out$s[i],'_lambda',out$lambda[k],'.score'), col.names=F, row.names=F, quote=F)
  }
}

#####
# Calculate the mean and sd of scores for each population specified in pop_scale
#####

# Calculate polygenic scores for reference individuals
# Do this using PLINK not lassosum as we want to use MAF to impute missing genotypes in target samples, which lassosum does not do.
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

for(i in 1:length(out$s)){
  for(k in 1:length(out$lambda)){
    system(paste0(opt$plink, ' --bfile ',opt$ref_plink_gw,' --score ',opt$output,'.s',out$s[i],'_lambda',out$lambda[k],'.score 2 4 6 sum --out ',opt$output,'.s',out$s[i],'_lambda',out$lambda[k],' --memory ',floor(opt$memory*0.7)))
  }
}

fam<-fread(paste0(opt$ref_plink_gw,'.fam'))
scores<-fam[,1:2]
names(scores)<-c('FID','IID')

for(i in 1:length(out$s)){
  for(k in 1:length(out$lambda)){
    profile<-fread(paste0(opt$output,'.s',out$s[i],'_lambda',out$lambda[k],'.profile'))
    scores<-cbind(scores, profile$SCORESUM)
    names(scores)[grepl('V2',names(scores))]<-paste0('SCORE_s',out$s[i],'_lambda',out$lambda[k])
  }
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

#####
# Calculate the mean and sd of scores for each population specified in pop_scale
#####
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

for(i in 1:length(out$s)){
  for(k in 1:length(out$lambda)){
    system(paste0('rm ',opt$output,'.s',out$s[i],'_lambda',out$lambda[k],'.profile'))
    system(paste0('rm ',opt$output,'.s',out$s[i],'_lambda',out$lambda[k],'.nosex'))
    system(paste0('rm ',opt$output,'.s',out$s[i],'_lambda',out$lambda[k],'.log'))
  }
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
