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
make_option("--sumstats", action="store", default=NA, type='character',
		help="GWAS summary statistics in LDSC format [required]"),
make_option("--gcta", action="store", default=NA, type='character',
		help="Path to GCTA binary [required]"),
make_option("--ldsc", action="store", default=NA, type='character',
		help="Path to LD-score regression binary [required]"),
make_option("--ldsc_ref", action="store", default=NA, type='character',
		help="Path to LD-score regression reference data 'eur_w_ld_chr' [required]"),
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
# Estimate the SNP-heritability
#####

system(paste0(opt$ldsc,' --h2 ',opt$sumstats,' --ref-ld-chr ',opt$ldsc_ref,'/ --w-ld-chr ',opt$ldsc_ref,'/ --out ', opt$output_dir,'ldsc_snp_h2_temp'))

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

GWAS_clean<-GWAS_clean[,c('V2','A1','A2','Z','P','N')]
names(GWAS_clean)<-c('SNP','A1','A2','Z','P','N')

nsnp<-dim(GWAS_clean)[1]

###
# Change to COJO format
###

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
GWAS_clean<-GWAS_clean[,c('SNP','A1','A2','Z','P','N','MAF')]

# Transform Z score to beta and se using formula from https://www.ncbi.nlm.nih.gov/pubmed/27019110
# Note, we could use full sumstats rather than munged which would contain more accurate beta and se.
GWAS_clean$beta<-GWAS_clean$Z/sqrt((2*GWAS_clean$MAF)*(1-GWAS_clean$MAF)*(GWAS_clean$N+sqrt(abs(GWAS_clean$Z))))
GWAS_clean$se<-abs(GWAS_clean$beta)/abs(GWAS_clean$Z)

GWAS_clean<-GWAS_clean[,c('SNP','A1','A2','MAF','beta','se','P','N'),with=F]
names(GWAS_clean)<-c('SNP','A1','A2','freq','b','se','p','N')

fwrite(GWAS_clean, paste0(opt$output_dir,'GWAS_sumstats_COJO.txt'), sep=' ')

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After harmonisation with the reference,',dim(GWAS_clean)[1],'variants remain.\n')
sink()

#####
# Run GCTA SBLUP
#####

system(paste0(opt$gcta,' --bfile ',opt$ref_plink,' --keep ',opt$ref_keep,' --cojo-file ',opt$output_dir,'GWAS_sumstats_COJO.txt --cojo-sblup ',nsnp*(1/ldsc_h2-1),' --cojo-wind 1000 --thread-num ',opt$n_cores,' --out ',opt$output_dir,'GWAS_sumstats_SBLUP'))

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
system(paste0('rm ',opt$output_dir,'GWAS_sumstats_COJO.txt'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
