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
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Number of cores for parallel computing [optional]"),
make_option("--sumstats", action="store", default=NA, type='character',
		help="GWAS summary statistics in LDSC format [required]"),
make_option("--nps", action="store", default=NA, type='character',
		help="Path to NPS directory [required]"),
make_option("--prune_hla", action="store", default=T, type='logical',
		help="Retain only top assocaited variant in HLA region [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

sink(file = paste(opt$output,'.log',sep=''), append = F)
cat(
'#################################################################
# polygenic_score_file_creator_NPS.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

#######
# Format the summary summary statistics
#######

# Read in sumstats
GWAS<-fread(cmd=paste0('zcat ',opt$sumstats))
GWAS<-GWAS[complete.cases(GWAS),]
GWAS<-GWAS[order(GWAS$CHR, GWAS$BP),]

# Convert OR to BETA
if(!('BETA' %in% names(GWAS))){
  GWAS$BETA<-log(GWAS$OR)
}

# Rename allele frequency column
if(sum(names(GWAS) == 'FREQ') != 1){
  GWAS$FREQ<-GWAS$REF.FREQ
}

# Flip FREQ to correspond to the reference allele
GWAS$FREQ<- 1-GWAS$FREQ

write.table(GWAS$SNP, paste0(opt$output_dir,'/GWAS_sumstats_snps.txt'), col.names=F, row.names=F, quote=F)

#######
# Create dosage format reference
#######

for(i in 1:22){
# Extract EUR individuals remove variants
system(paste0(opt$plink,' --bfile ',opt$ref_plink_chr,i,' --extract ',opt$output_dir,'/GWAS_sumstats_snps.txt --make-bed --keep ',opt$ref_keep,' --out ',opt$output_dir,'/ref_chr',i))
system(paste0(opt$plink,' --bfile ',opt$ref_plink_chr,i,' --extract ',opt$output_dir,'/GWAS_sumstats_snps.txt --recode vcf --maf 0.01 --keep ',opt$ref_keep,' --geno 0 --out ',opt$output_dir,'/ref_chr',i))
  
# Convert to dosage
system(paste0(opt$qctool,' -g ',opt$output_dir,'/ref_chr',i,'.vcf -og ',opt$output_dir,'/ref_chr',i,'.dosage'))

system(paste0('rm ',opt$output_dir,'/ref_chr',i,'.vcf'))
system(paste0('rm ',opt$output_dir,'/ref_chr',i,'.nosex'))
system(paste0('rm ',opt$output_dir,'/ref_chr',i,'.log'))

system(paste0('/users/k1806347/brc_scratch/Software/pigz ',opt$output_dir,'/ref_chr',i,'.dosage'))

system(paste0('mv ',opt$output_dir,'/ref_chr',i,'.dosage.gz ',opt$output_dir,'/chrom',i,'.temp.dosage.gz'))
}

#######
# Harmonise alleles with the reference
#######

ref_bim<-NULL
for(i in 1:22){
  ref_bim<-rbind(ref_bim, fread(paste0(opt$output_dir,'/ref_chr',i,'.bim')))
}

GWAS_match<-merge(GWAS, ref_bim[,c('V2','V5','V6'),with=F], by.x=c('SNP','A1','A2'), by.y=c('V2','V5','V6'))
GWAS_switch<-merge(GWAS, ref_bim[,c('V2','V5','V6'),with=F], by.x=c('SNP','A1','A2'), by.y=c('V2','V6','V5'))
GWAS_switch$allele_tmp<-GWAS_switch$A2
GWAS_switch$A2<-GWAS_switch$A1
GWAS_switch$A1<-GWAS_switch$allele_tmp
GWAS_switch$allele_tmp<-NULL
GWAS_switch$BETA<--GWAS_switch$BETA
GWAS_switch$FREQ<-1-GWAS_switch$FREQ
GWAS<-rbind(GWAS_match, GWAS_switch)

GWAS[GWAS$SNP == 'rs12027409',]

# Subset column and rename
GWAS<-GWAS[,c('CHR','BP','A2','A1','FREQ','P','BETA'), with=F]
names(GWAS)<-c('chr','pos','ref','alt','reffreq','pval','effalt')

GWAS$chr<-paste0('chr',GWAS$chr)

fwrite(GWAS, paste0(opt$output_dir,'/GWAS_sumstats_tmp.txt'), na='NA', sep='\t', quote=F)

#######
# Standardised genotypes
#######

setwd(opt$nps)
system(paste0('./run_all_chroms.sh sge/nps_stdgt.job ',opt$output_dir,' temp'))

#######
# Configure NPS run
#######

setwd(opt$nps)
system(paste0('Rscript npsR/nps_init.R --gwas ',opt$output_dir,'/GWAS_sumstats_tmp.txt --train-fam ',opt$output_dir,'/ref_chr22.fam --train-dir ',opt$output_dir,' --train-dataset temp --out ',opt$output_dir,'/config'))

# NOTE. this step does not work because there is no phenotypic data in the reference file.
# This distringuishes this method from others, as a score file cannot be generated without using a sample with phenotype data.

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
