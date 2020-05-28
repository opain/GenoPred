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
make_option("--gctb", action="store", default=NA, type='character',
    help="Path to GCTB binary [required]"),
make_option("--munged", action="store", default=T, type='logical',
    help="Logical indicating whether the GWAS summary statistics are in munged format [required]"),
make_option("--impute_N", action="store", default=T, type='logical',
    help="Logical indicating whether per variant N should imputed based on SE. Ony possible when munged == F [optional]"),
make_option("--auto_trunc", action="store", default=F, type='logical',
    help="Logical indicating whether optimal sample size truncation should be used. Default is 0.1*median. [required]"),
make_option("--P_max", action="store", default=NA, type='numeric',
    help="P-value threshold for filter variants [optional]"),
make_option("--ld_matrix_chr", action="store", default=NA, type='character',
		help="Path to per chromosome shrunk sparse LD matrix from GCTB [required]")
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
# Read in sumstats and insert p-values
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in GWAS and harmonising with reference.\n')
sink()

GWAS<-fread(cmd=paste0('zcat ',opt$sumstats))
GWAS<-GWAS[complete.cases(GWAS),]

# If unmunged convert OR to BETA, and do some basic QC
if(opt$munged == F){
  if(!('BETA' %in% names(GWAS))){
    GWAS$BETA<-log(GWAS$OR)
  }
  if(('FREQ' %in% names(GWAS))){
    GWAS<-GWAS[GWAS$FREQ>0.01 & GWAS$FREQ<0.99,]
  }
  if(('INFO' %in% names(GWAS))){
    GWAS<-GWAS[GWAS$INFO>0.6,]
  }
}

if(opt$munged == F){
  if(!('SE' %in% names(GWAS))){
    GWAS$Z<-abs(qnorm(GWAS$P/2))
    GWAS$SE<-abs(GWAS$BETA/GWAS$Z)
    GWAS<-GWAS[!is.na(GWAS$SE) & GWAS$SE != 0,]
  }
}

if(opt$munged == T){
  # Calculate P values (first change 0 z-scores to a small non-zero number, as otherwise calculating the beta and se leads to na)
  GWAS$Z[GWAS$Z == 0]<-1e-10
  GWAS$P<-2*pnorm(-abs(GWAS$Z))
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('GWAS contains',dim(GWAS)[1],'variants.\n')
sink()

if(!is.na(opt$P_max)){
  GWAS<-GWAS[GWAS$P < opt$P_max,]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('GWAS contains',dim(GWAS)[1],'variants after removing variants with P >=',opt$P_max,'\n')
  sink()
}

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

if(opt$munged == T){
  GWAS_clean<-GWAS_clean[,c('V2','A1','A2','Z','P','N')]
  names(GWAS_clean)<-c('SNP','A1','A2','Z','P','N')
} else {
  if(!('FREQ' %in% names(GWAS_clean))){
    GWAS_clean<-GWAS_clean[,c('V2','A1','A2','BETA','SE','P','N')]
    names(GWAS_clean)<-c('SNP','A1','A2','BETA','SE','P','N')
  } else {
    GWAS_clean<-GWAS_clean[,c('V2','A1','A2','BETA','SE','P','N','FREQ')]
    names(GWAS_clean)<-c('SNP','A1','A2','BETA','SE','P','N','MAF')
  }
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('After merging with the reference there are',dim(GWAS_clean)[1],'variants.\n')
sink()

###
# Change to COJO format
###

if(!('MAF' %in% names(GWAS_clean))){
  if(opt$munged == F){
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Using reference MAF estimates.\n')
    sink()
  }
    
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
  if(opt$munged == T){
    GWAS_clean<-GWAS_clean[,c('SNP','A1','A2','Z','P','N','MAF')]
  } else {
    GWAS_clean<-GWAS_clean[,c('SNP','A1','A2','BETA','SE','P','N','MAF')]
  }
}

# Truncate the sample size distribution across variants
if(length(unique(GWAS_clean$N)) == 1){
	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat('GWAS sumstats do not include variant specific sample size info.\n')
	sink()
	
	per_var_N<-F

} else {

  per_var_N<-T
  
	bitmap(paste(opt$output,'N_dist_preTrunc.png',sep=''), units='px', res=300, width=3000, height=5000)
	hist(GWAS_clean$N)
	dev.off()

	med_N<-median(GWAS_clean$N)

	if(opt$auto_trunc == F){
	  GWAS_clean<-GWAS_clean[GWAS_clean$N > (med_N-(med_N*0.10)) & GWAS_clean$N < (med_N+(med_N*0.10)),]
	  
	} else {
	  N_sd<-NULL
	  sd_dev<-seq(0,0.3,0.005)
	  for(i in 1:length(sd_dev)){
	    N_sd<-rbind(N_sd, data.frame(i=sd_dev[i],
	                                 sd_i=sd(scale(GWAS_clean$N)[GWAS_clean$N > (med_N-(med_N*sd_dev[i])) & GWAS_clean$N < (med_N+(med_N*sd_dev[i]))]),
	                                 sd_diff=NA))
        if(i!=1){
  	      N_sd$sd_diff[i]<-N_sd$sd_i[i] - N_sd$sd_i[i-1]
  	    }
	  }
	  
	  opt_i<-N_sd$i[which(N_sd$sd_diff==max(N_sd$sd_diff,na.rm=T))-1]
	  GWAS_clean<-GWAS_clean[GWAS_clean$N > (med_N-(med_N*opt_i)) & GWAS_clean$N < (med_N+(med_N*opt_i)),]
	
	  sink(file = paste(opt$output,'.log',sep=''), append = T)
	  cat('Variants removed is N if ',round(opt_i*100,2),'% from the median N.\n', sep='')
	  sink()
	}

	bitmap(paste(opt$output,'N_dist_postTrunc.png',sep=''), units='px', res=300, width=3000, height=5000)
	hist(GWAS_clean$N)
	dev.off()

	sink(file = paste(opt$output,'.log',sep=''), append = T)
	cat(dim(GWAS_clean)[1],'variants in GWAS sumstats after truncating sample size distribution.\n')
	sink()
}

if(opt$munged == T){
  # Transform Z score to beta and se using formula from https://www.ncbi.nlm.nih.gov/pubmed/27019110
  # Note, we could use full sumstats rather than munged which would contain more accurate beta and se.
  GWAS_clean$beta<-GWAS_clean$Z/sqrt((2*GWAS_clean$MAF)*(1-GWAS_clean$MAF)*(GWAS_clean$N+sqrt(abs(GWAS_clean$Z))))
  GWAS_clean$se<-abs(GWAS_clean$beta)/abs(GWAS_clean$Z)
  
  GWAS_clean<-GWAS_clean[,c('SNP','A1','A2','MAF','beta','se','P','N'),with=F]
} else {
  GWAS_clean<-GWAS_clean[,c('SNP','A1','A2','MAF','BETA','SE','P','N'),with=F]
}

GWAS_clean<-GWAS_clean[GWAS_clean$MAF != 0,]
GWAS_clean<-GWAS_clean[GWAS_clean$MAF != 1,]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat(dim(GWAS_clean)[1],'variants in GWAS sumstats after removing invariant variation.\n')
sink()

names(GWAS_clean)<-c('SNP','A1','A2','freq','b','se','p','N')

fwrite(GWAS_clean, paste0(opt$output_dir,'GWAS_sumstats_COJO.txt'), sep=' ', na = "NA", quote=F)

#####
# Run GCTB SBayesR
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Running SBayesR analysis...')
sink()

is.odd <- function(x){x %% 2 != 0}
CHROMS<-c(1:22, rep(NA,30)) 
CHROMS_mat<-matrix(CHROMS,nrow=opt$n_cores, ncol=ceiling(22/opt$n_cores)) 
for(i in which(is.odd(1:dim(CHROMS_mat)[2]))){CHROMS_mat[,i]<-rev(CHROMS_mat[,i])} 
CHROMS<-as.numeric(CHROMS_mat) 
CHROMS<-CHROMS[!is.na(CHROMS)] 
print(CHROMS)

error<-foreach(i=CHROMS, .combine=rbind) %dopar% {
  if(per_var_N == F & opt$impute_N == T){
    log<-system(paste0(opt$gctb,' --sbayes R --ldm ',opt$ld_matrix_chr,i,'.ldm.sparse --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1 --gwas-summary ',opt$output_dir,'GWAS_sumstats_COJO.txt --chain-length 10000 --exclude-mhc --burn-in 2000 --impute-n --out-freq 1000 --out ',opt$output_dir,'GWAS_sumstats_SBayesR.chr',i), intern=T)
  } else {
    log<-system(paste0(opt$gctb,' --sbayes R --ldm ',opt$ld_matrix_chr,i,'.ldm.sparse --pi 0.95,0.02,0.02,0.01 --gamma 0.0,0.01,0.1,1 --gwas-summary ',opt$output_dir,'GWAS_sumstats_COJO.txt --chain-length 10000 --exclude-mhc --burn-in 2000 --out-freq 1000 --out ',opt$output_dir,'GWAS_sumstats_SBayesR.chr',i), intern=T)
  }
  
	if(sum(grepl("MCMC cycles completed", log) == T) == 1 & sum(grepl("Analysis finished", log) == T) == 1){
	  error_log<-data.frame(chr=i, Log='Analysis converged')
	}
	
	if(sum(grepl("MCMC cycles completed", log) == T) == 0 & sum(grepl("Analysis finished", log) == T) == 1){
	  error_log<-data.frame(chr=i, Log='Analysis did not converge')
	}
	
	if(sum(grepl("MCMC cycles completed", log) == T) == 0 & sum(grepl("Analysis finished", log) == T) == 0){
	  error_log<-data.frame(chr=i, Log='Error')
	}
	
	error_log
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

if(sum(grepl('Error', error$Log) == T) > 1){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('An error occurred for',sum(grepl('Error', error$Log) == T),'chromosomes. Retry requesting more memory or run interactively to debug.\n')
    print(error)
  sink()
  system(paste0('rm ',opt$output_dir,'GWAS_sumstats_*'))
  q()
  n
}

# Check whether analysis completed for all chromosomes
comp_list<-list.files(path=opt$output_dir, pattern='GWAS_sumstats_SBayesR.chr.*snpRes')
incomp<-c(1:22)[!(paste0('GWAS_sumstats_SBayesR.chr',1:22,'.snpRes') %in% comp_list)]
comp<-c(1:22)[(paste0('GWAS_sumstats_SBayesR.chr',1:22,'.snpRes') %in% comp_list)]

# Combine per chromosome snpRes files
snpRes<-NULL
for(i in comp){
  snpRes<-rbind(snpRes, fread(paste0(opt$output_dir,'GWAS_sumstats_SBayesR.chr',i,'.snpRes')))
}

snpRes<-snpRes[,c('Chrom','Name','Position','A1','A2','A1Effect')]

write.table(snpRes, paste0(opt$output_dir,'GWAS_sumstats_SBayesR.GW.snpRes'), col.names=F, row.names=F, quote=F)

# Combine per chromosome parRes files
parRes_mcmc<-list()
for(i in comp){
	parRes_mcmc[[i]]<-fread(paste0(opt$output_dir,'GWAS_sumstats_SBayesR.chr',i,'.mcmcsamples.Par'))
}

parRes<-NULL
for(par in names(parRes_mcmc[[i]])){
	parRes_mcmc_par<-NULL
	for(i in comp){
		parRes_mcmc_par<-cbind(parRes_mcmc_par,parRes_mcmc[[i]][[par]])
	}
	
	parRes_mcmc_par_sum<-rowSums(parRes_mcmc_par)
	
	parRes_par<-data.frame(	Par=par,
							Mean=mean(parRes_mcmc_par_sum),
							SD=sd(parRes_mcmc_par_sum))
							
	parRes<-rbind(parRes,parRes_par)
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('SNP-heritability estimate is ',parRes[parRes$Par == 'hsq', names(parRes) == 'Mean']," (SD=",parRes[parRes$Par == 'hsq', names(parRes) == 'SD'],").\n", sep='')
sink()

write.table(parRes, paste0(opt$output_dir,'GWAS_sumstats_SBayesR.GW.parRes'), col.names=T, row.names=F, quote=F)

system(paste0('rm ',opt$output_dir,'GWAS_sumstats_SBayesR.chr*')) 
system(paste0('rm ',opt$output_dir,'GWAS_sumstats_COJO.txt')) 

sink(file = paste(opt$output,'.log',sep=''), append = T)
if(length(incomp) >=1){
	cat('No results for chromosomes',paste(incomp, collapse=','),'\n')
}
cat('Logs:\n')
print(error)
cat('\n')
sink()

####
# Calculate mean and sd of polygenic scores at each threshold
####

# Calculate polygenic scores for reference individuals
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

system(paste0(opt$plink, ' --bfile ',opt$ref_plink,' --score ',opt$output_dir,'GWAS_sumstats_SBayesR.GW.snpRes 2 4 6 sum --out ',opt$output_dir,'ref.profiles --memory ',floor(opt$memory*0.7)))

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

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Clean up temporary files
###

system(paste0('rm ',opt$output_dir,'ref.profiles.*'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
