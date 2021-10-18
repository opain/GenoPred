#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--ref_plink", action="store", default=NA, type='character',
    help="Path to per chromosome reference PLINK files [required]"),
make_option("--ldpred2_ref_dir", action="store", default=NA, type='character',
    help="Path to directory containing LDPred2 reference data [required]"),
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
make_option("--sd_check", action="store", default=T, type='logical',
    help="Set to T to compare GWAS SD to reference SD [optional]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--binary", action="store", default=F, type='logical',
    help="Specify T if GWAS phenotyp is binary [optional]"),
make_option("--sumstats", action="store", default=NA, type='character',
		help="GWAS summary statistics [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)
library(bigsnpr)

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
# polygenic_score_file_creator_LDPred2.R V1.0
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

# This script uses instructions from https://privefl.github.io/bigsnpr/articles/LDpred2.html

#####
# Read in reference data
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in reference.\n')
sink()

# Attach the "bigSNP" object in R session
if(file.exists(paste0(opt$output_dir,'ref.bk'))){
  system(paste0('rm ', opt$output_dir,'ref.bk'))
}

snp_readBed(paste0(opt$ref_plink,'.bed'), backingfile = paste0(opt$output_dir,'ref'))
ref <- snp_attach(paste0(opt$output_dir,"ref.rds"))

ref_keep<-read.table(opt$ref_keep, header=F, stringsAsFactors=F)
ind_keep<-which(ref$fam$family.ID %in% ref_keep$V1)

G   <- ref$genotypes
CHR <- ref$map$chromosome
POS <- ref$map$physical.pos
y   <- ref$fam$affection - 1
NCORES <- opt$n_cores

#####
# Read in sumstats
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in GWAS.\n')
sink()

sumstats <- bigreadr::fread2(opt$sumstats)
str(sumstats)

sumstats<-sumstats[complete.cases(sumstats),]
sumstats<-sumstats[sumstats$SE != 0,]

# Extract subset if testing
if(!is.na(opt$test)){
  if(single_chr_test == F){
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Testing mode enabled. Extracted ',opt$test,' variants per chromsome.\n', sep='')
    sink()
    
    GWAS_test<-NULL
    for(i in 1:22){
      GWAS_tmp<-sumstats[sumstats$CHR == i,]
      GWAS_tmp<-GWAS_tmp[order(GWAS_tmp$BP),]
      GWAS_tmp<-GWAS_tmp[1:opt$test,]
      GWAS_test<-rbind(GWAS_test,GWAS_tmp)
    }
    
    sumstats<-GWAS_test
    sumstats<-sumstats[complete.cases(sumstats),]
    rm(GWAS_test)
    print(table(sumstats$CHR))
    
  } else {
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Testing mode enabled. Extracted chromosome ',opt$test,' variants per chromsome.\n', sep='')
    sink()
    
    sumstats<-sumstats[sumstats$CHR == CHROMS,]
    print(table(sumstats$CHR))
  }
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Sumstats contains',dim(sumstats)[1],'variants.\n')
sink()

# Convert OR into BETA
if(sum(names(sumstats) == 'OR') == 1){
  sumstats$BETA<-log(sumstats$OR)
}

# Format sumstats for as in LDPred2 tutorial
if(sum(names(sumstats) == 'Ncas') == 1 & sum(names(sumstats) == 'Ncon') == 1){
  sumstats$n_eff <- 4 / (1 / sumstats$Ncas + 1 / sumstats$Ncon)
  sumstats$Ncas <- sumstats$Ncon <- NULL
  sumstats<-sumstats[,c('CHR','SNP','BP','A1','A2','BETA','SE','n_eff','P')]
  names(sumstats)<-c('chr','rsid','pos','a1','a0','beta','beta_se','n_eff','p')
} else {
  sumstats<-sumstats[,c('CHR','SNP','BP','A1','A2','BETA','SE','N','P')]
  names(sumstats)<-c('chr','rsid','pos','a1','a0','beta','beta_se','n_eff','p')
}

if(!is.na(opt$test)){
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  test_start.time <- Sys.time()
  cat('Test started at',as.character(test_start.time),'\n')
  sink()
}

# Harmonise with the reference
map <- ref$map[-(2:3)]
names(map) <- c("chr", "pos", "a1", "a0")
info_snp <- snp_match(sumstats, map)

#####
# Perform additional suggested QC for LDPred2
#####

if(opt$sd_check == T){
  # Read in reference genotype SD
  sd<-readRDS(paste0(opt$ldpred2_ref_dir,'/sd.rds'))
  
  # Remove SDss<0.5???SDval or SDss>0.1+SDval or SDss<0.1 or SDval<0.05
  sd_val <- sd[info_snp$`_NUM_ID_`]

  if(opt$binary == F){
    sd_y_est = median(sd_val * info_snp$beta_se * sqrt(info_snp$n_eff))
    sd_ss = with(info_snp, sd_y_est / sqrt(n_eff * beta_se^2))
  } else {
    sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))
  }
  
  is_bad <-sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
  
  library(ggplot2)
  bitmap(paste0(opt$output_dir,'/LDPred2_sd_qc.png'), res=300, unit='px',height=2000, width=2000)
  print(qplot(sd_val, sd_ss, color = is_bad) +
    theme_bigstatsr() +
    coord_equal() +
    scale_color_viridis_d(direction = -1) +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Standard deviations in the validation set",
         y = "Standard deviations derived from the summary statistics",
         color = "Removed?"))
  dev.off()
  
  # If more than half the variants have the wrong SD then the N is probably inaccurate
  # Recompute N based on BETA and SE
  if(sum(is_bad) > (length(is_bad)*0.5)){
    n_eff_imp<-sd_val^2/info_snp$beta_se^2
    
    n_eff_imp<-((2/sd_val)^2)/(info_snp$beta_se^2)
    info_snp$n_eff<-median(n_eff_imp)
    
    sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2))
    is_bad <-sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
    
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    test_start.time <- Sys.time()
    cat('More than half the variants had a discordant SD.\n')
    cat('Median imputed N estimate: ', median(n_eff_imp),'.\n',sep='')
    sink()
    
    library(ggplot2)
    bitmap(paste0(opt$output_dir,'/LDPred2_sd_qc_impN.png'), res=300, unit='px',height=2000, width=2000)
    print(qplot(sd_val, sd_ss, color = is_bad) +
      theme_bigstatsr() +
      coord_equal() +
      scale_color_viridis_d(direction = -1) +
      geom_abline(linetype = 2, color = "red") +
      labs(x = "Standard deviations in the validation set",
           y = "Standard deviations derived from the summary statistics",
           color = "Removed?"))
    dev.off()
    
  }
  
  sumstats<-info_snp[!is_bad, ]
  
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Sumstats contains',dim(sumstats)[1],'after additional genotype SD check.\n')
  sink()
}

#######
# Estimate heritability
#######

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Reading in reference LD scores.\n')
sink()

# Read in LD score estimates
ld_score<-readRDS(paste0(opt$ldpred2_ref_dir,'/map.rds'))
ld_score<-ld_score[info_snp$`_NUM_ID_`,]
ld_score<-ld_score[!is_bad, ]

ldsc <- with(sumstats, snp_ldsc(ld_score$ld, length(ld_score$ld), chi2 = (beta / beta_se)^2, sample_size = n_eff, blocks = NULL))

h2_est <- ldsc[["h2"]]

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Estimated SNP-based heritability =',h2_est,'\n')
sink()

if(!is.na(opt$test) & h2_est < 0.05){
  h2_est<-0.05
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Creating genome-wide sparse matrix.\n')
sink()

# Create genome-wide sparse LD matrix
for (chr in CHROMS) {
  print(chr)
  
  ## indices in 'sumstats'
  ind.chr <- which(sumstats$chr == chr)
  ## indices in 'map'
  ind.chr2 <- sumstats$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map$chr == chr))
  
  corr0 <- readRDS(paste0(opt$ldpred2_ref_dir,'/LD_chr', chr, ".rds"))[ind.chr3, ind.chr3]

  if (chr == CHROMS[1]) {
    corr <- as_SFBM(corr0, paste0(opt$output_dir,'/LD_GW_sparse'), compact = TRUE)
  } else {
    corr$add_columns(corr0, nrow(corr))
  }
}

#####
# Run LDPred2
#####

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Running LDPred models.\n')
sink()

# LDpred2-inf
beta_inf <- snp_ldpred2_inf(corr, sumstats, h2_est)

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Infintesimal model complete at',as.character(Sys.time()),'\n')
sink()

# LDpred2-grid
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

beta_grid <- snp_ldpred2_grid(corr, sumstats, params, ncores = NCORES)

beta_grid_nosp<-data.table(beta_grid[,params$sparse == F])
names(beta_grid_nosp)<-gsub('-','.',paste0(params$p[params$sparse == F],'_',params$h2[params$sparse == F],'_nosparse'))

beta_grid_sp<-data.table(beta_grid[,params$sparse == T])
names(beta_grid_sp)<-gsub('-','.',paste0(params$p[params$sparse == T],'_',params$h2[params$sparse == T],'_sparse'))

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Grid model complete at',as.character(Sys.time()),'\n')
sink()

# LDpred2-auto
multi_auto <- snp_ldpred2_auto(corr, sumstats, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = NCORES)

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

pred_auto <- big_prodMat(G, beta_auto, ind.row = ind_keep,
                         ind.col = sumstats[["_NUM_ID_"]])

sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_beta_auto <- rowMeans(beta_auto[, keep])

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Auto model complete at',as.character(Sys.time()),'\n')
sink()

# compute predictions for test set
betas <- data.table(SNP=sumstats$rsid, A1=sumstats$a1, A2=sumstats$a0, beta_inf, beta_grid_nosp, beta_grid_sp, beta_auto = final_beta_auto)
names(betas)[-1:-3]<-paste0('SCORE_',names(betas)[-1:-3])

rem<-NULL
for(i in 4:length(names(betas))){
  if(is.infinite(sum(betas[[names(betas)[i]]])) | is.na(sum(betas[[names(betas)[i]]]))){
    sink(file = paste(opt$output,'.log',sep=''), append = T)
    cat('Skipping',names(betas)[i],'due to presence of non-finite values.\n')
    sink()
    rem<-c(rem,i)
  } else {
    fwrite(betas[,c('SNP','A1','A2',names(betas)[i]), with=F], paste0(opt$output,'.',gsub('SCORE_','',names(betas)[i]),'.SCORE'), quote=F, sep='\t', col.names=F, na='NA')
  }
}

if(is.null(rem) == F){
  betas<-betas[,-rem, with=F]
}

if(!is.na(opt$test)){
  end.time <- Sys.time()
  time.taken <- end.time - test_start.time
  sink(file = paste(opt$output,'.log',sep=''), append = T)
  cat('Test run finished at',as.character(end.time),'\n')
  cat('Test duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
  sink()
  system(paste0('rm ',opt$output_dir,'*.SCORE'))
  system(paste0('rm ',opt$output_dir,'LD_GW_sparse.sbk'))
  system(paste0('rm ',opt$output_dir,'ref.rds'))
  system(paste0('rm ',opt$output_dir,'ref.bk'))
  q()
}

####
# Calculate mean and sd of polygenic scores at each threshold
####

# Calculate polygenic scores for reference individuals
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Calculating polygenic scores in reference...')
sink()

param<-gsub('SCORE_','',names(betas)[-1:-3])

for(i in param){
  system(paste0(opt$plink, ' --bfile ',opt$ref_plink,' --score ',opt$output,'.',i,'.SCORE 1 2 4 sum --out ',opt$output_dir,'ref.profiles.',i,' --memory ',floor(opt$memory*0.7)))
}

fam<-fread(paste0(opt$ref_plink,'.fam'))
scores<-fam[,1:2]
names(scores)<-c('FID','IID')

for(i in param){
  SCORE_temp<-fread(paste0(opt$output_dir,'ref.profiles.',i,'.profile'))
  scores<-cbind(scores, SCORE_temp[, 6])
  names(scores)[grepl('SCORESUM',names(scores))]<-paste0('SCORE_',i)
}

# Calculate the mean and sd of scores for each population specified in pop_scale
pop_keep_files<-read.table(opt$ref_pop_scale, header=F, stringsAsFactors=F)

for(k in 1:dim(pop_keep_files)[1]){
  pop<-pop_keep_files$V1[k]
  keep<-fread(pop_keep_files$V2[k], header=F)
  scores_keep<-scores[(scores$FID %in% keep$V1),]
  
  ref_scale<-data.frame(	Param=names(scores_keep[,-1:-2]),
                         Mean=sapply(scores_keep[,-1:-2], function(x) mean(x)),
                         SD=sapply(scores_keep[,-1:-2], function(x) sd(x)))
  
  fwrite(ref_scale, paste0(opt$output,'.',pop,'.scale'), sep=' ')
}

sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Done!\n')
sink()

###
# Clean up temporary files
###

system(paste0('rm ',opt$output_dir,'ref.profiles.*'))
system(paste0('rm ',opt$output_dir,'LD_GW_sparse.sbk'))
system(paste0('rm ',opt$output_dir,'ref.rds'))
system(paste0('rm ',opt$output_dir,'ref.bk'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
