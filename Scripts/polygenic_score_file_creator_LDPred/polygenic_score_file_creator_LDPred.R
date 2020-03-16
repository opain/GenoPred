#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--ref_plink", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK files [required]"),
make_option("--ref_keep", action="store", default=NA, type='character',
		help="Keep file to subset individuals in reference for clumping [required]"),
make_option("--ref_freq_chr", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK .frq files [required]"),
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
make_option("--sumstats", action="store", default=NA, type='character',
		help="GWAS summary statistics in LDSC format [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(foreach)
library(doMC)
registerDoMC(opt$n_cores)

tmp<-sub('.*/','',opt$output)
opt$output_dir<-sub(paste0(tmp,'*.'),'',opt$output)
system(paste0('mkdir -p ',opt$output_dir))

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

# Check GWAS_N
GWAS_N<-mean(GWAS$N)

# Calculate P values
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
bim<-fread(paste0(opt$ref_plink,'.bim'))

bim$IUPAC[bim$V5 == 'A' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='A']<-'W'
bim$IUPAC[bim$V5 == 'C' & bim$V6 =='G' | bim$V5 == 'G' & bim$V6 =='C']<-'S'
bim$IUPAC[bim$V5 == 'A' & bim$V6 =='G' | bim$V5 == 'G' & bim$V6 =='A']<-'R'
bim$IUPAC[bim$V5 == 'C' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='C']<-'Y'
bim$IUPAC[bim$V5 == 'G' & bim$V6 =='T' | bim$V5 == 'T' & bim$V6 =='G']<-'K'
bim$IUPAC[bim$V5 == 'A' & bim$V6 =='C' | bim$V5 == 'C' & bim$V6 =='A']<-'M'

bim_GWAS<-merge(bim,GWAS, by.x='V2', by.y='SNP')
GWAS_clean<-bim_GWAS[bim_GWAS$IUPAC.x == bim_GWAS$IUPAC.y,]

GWAS_clean<-GWAS_clean[,c('V1','V2','V4','A1','A2','Z','P','N')]
names(GWAS_clean)<-c('CHR','SNP','BP','A1','A2','Z','P','N')

# Insert frq of each variant based on reference data
freq<-NULL
for(i in 1:22){
  freq_tmp<-fread(paste0(opt$ref_freq_chr,i,'.frq'))
  freq<-rbind(freq, freq_tmp)
}

GWAS_clean_frq_match<-merge(GWAS_clean, freq, by=c('SNP','A1','A2'))
GWAS_clean_frq_switch<-merge(GWAS_clean, freq, by.x=c('SNP','A1','A2'), by.y=c('SNP','A2','A1'))
GWAS_clean_frq_switch$MAF<-1-GWAS_clean_frq_switch$MAF
GWAS_clean<-rbind(GWAS_clean_frq_match, GWAS_clean_frq_switch)
GWAS_clean<-GWAS_clean[,c('CHR.x','SNP','BP','A1','A2','Z','P','N','MAF')]
names(GWAS_clean)<-c('CHR','SNP','BP','A1','A2','Z','P','N','MAF')

# Flip MAF to be reference allele freq
GWAS_clean$MAF<-1-GWAS_clean$MAF

# Remove invariant SNPs
GWAS_clean<-GWAS_clean[GWAS_clean$MAF != 0,]
GWAS_clean<-GWAS_clean[GWAS_clean$MAF != 1,]

# Transform Z score to beta using formula from https://www.ncbi.nlm.nih.gov/pubmed/27019110
GWAS_clean$beta<-GWAS_clean$Z/sqrt((2*GWAS_clean$MAF)*(1-GWAS_clean$MAF)*(GWAS_clean$N+sqrt(abs(GWAS_clean$Z))))

# Insert fake INFO score
GWAS_clean$INFO<-1

GWAS_clean<-GWAS_clean[,c('CHR','BP','A1','A2','MAF','INFO','SNP','P','beta'),with=F]
names(GWAS_clean)<-c('chr','pos','alt','ref','reffrq','info','rs','pval','effalt')

GWAS_clean$chr<-paste0('chr',GWAS_clean$chr)

fwrite(GWAS_clean, paste0(opt$output_dir,'GWAS_sumstats.txt'), sep=' ', na = "NA", quote=F)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After harmonisation with the reference,',dim(GWAS_clean)[1],'variants remain.\n')
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

#####
# Corrdinate the GWAS summary statistics with reference using LDPred
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Corrdinating GWAS summary statistics with reference using LDPred...')
sink()

log<-system(paste0(opt$ldpred,' coord --ssf-format STANDARD --N ',round(GWAS_N,0),' --ssf ',opt$output_dir,'GWAS_sumstats.txt --out ',opt$output_dir,'GWAS_sumstats.coord --gf ',opt$ref_plink_subset), intern=T)

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

system(paste0(opt$ldpred,' gibbs --cf ',opt$output_dir,'GWAS_sumstats.coord --ldr ',round(nsnp/3000,0),' --ldf ', opt$ref_plink_subset,' --out ',opt$output,'.weight --N ',round(GWAS_N,0)))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

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
