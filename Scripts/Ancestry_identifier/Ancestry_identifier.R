#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--target_plink_chr", action="store", default=NULL, type='character',
    help="Path to per chromosome target PLINK2 files [required]"),
make_option("--target_keep", action="store", default=NULL, type='character',
    help="Path to file listing individuals in the target sample to retain [optional]"),
make_option("--ref_plink_chr", action="store", default=NULL, type='character',
    help="Path to per chromosome reference PLINK2 files [required]"),
make_option("--ref_keep", action="store", default=NULL, type='character',
    help="Path to file listing individuals in the reference sample to retain [optional]"),
make_option("--maf", action="store", default=0.05, type='numeric',
    help="Minor allele frequency threshold [optional]"),
make_option("--geno", action="store", default=0.02, type='numeric',
    help="Variant missingness threshold [optional]"),
make_option("--hwe", action="store", default=1e-6, type='numeric',
    help="Hardy Weinberg p-value threshold. [optional]"),
make_option("--n_pcs", action="store", default=6, type='numeric',
		help="Number of PCs (min=4) [optional]"),
make_option("--plink2", action="store", default='plink2', type='character',
		help="Path PLINK software binary [optional]"),
make_option("--output", action="store", default=NULL, type='character',
		help="Path for output files [required]"),
make_option("--pop_data", action="store", default=NULL, type='character',
    help="Population data for the reference samples [required]"),
make_option("--model_method", action="store", default='glmnet', type='character',
    help="Method used for generate prediction model [optional]"),
make_option("--sd_rule", action="store", default=F, type='logical',
    help="Logical indicating whether the 3SD rule should be used to define ancestry, or the model-based approach [optional]"),
make_option("--prob_thresh", action="store", default=0.95, type='numeric',
    help="Indicates whether probability threshold should be used when defining ancestry [optional]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify test mode [optional]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
source('../functions/misc.R')
source_all('../functions')
library(data.table)
library(caret)
library(pROC)
library(verification)
library(ggplot2)
library(cowplot)

# Check required inputs
if(is.null(opt$target_plink_chr)){
  stop('--target_plink_chr must be specified.\n')
}
if(is.null(opt$ref_plink_chr)){
  stop('--ref_plink_chr must be specified.\n')
}
if(is.null(opt$output)){
  stop('--output must be specified.\n')
}
if(is.null(opt$pop_data)){
  stop('--pop_data must be specified.\n')
}

# Create output directory
opt$out_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$out_dir))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'Ancestry_identifier.R', start.time = start.time)


# If testing, change CHROMS to chr value, and lower ancestry probability threshold
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
  opt$prob_thresh <- 0.5
  log_add(log_file = log_file, message = 'Lowering prob_thresh parameter to 0.5 for testing.')
}

if(nrow(fread(paste0(opt$ref_plink_chr, CHROMS[1],'.psam'))) < 100){
  stop('opt$ref_plink_chr must contain at least 100 individuals.')
}

###########
# Extract target_keep
###########

if(!is.null(opt$target_keep)){
  plink_subset(keep = opt$target_keep, chr = CHROMS, plink2 = opt$plink2, pfile = opt$target_plink_chr, out = paste0(tmp_dir,'/target_subset.chr'))
  opt$target_plink_chr_subset<-paste0(tmp_dir,'/target_subset')
} else {
  opt$target_plink_chr_subset<-opt$target_plink_chr
}

###########
# Extract ref_keep
###########

if(!is.null(opt$ref_keep)){
  plink_subset(keep = opt$ref_keep, chr = CHROMS, plink2 = opt$plink2, pfile = opt$ref_plink_chr, out = paste0(tmp_dir,'/ref_subset.chr'))
  opt$ref_plink_chr_subset<-paste0(tmp_dir,'/ref_subset.chr')
} else {
  opt$ref_plink_chr_subset<-opt$ref_plink_chr
}

###########
# QC target
###########

# If target sample size is <100, only apply SNP missingness parameter
psam<-fread(paste0(opt$target_plink_chr_subset, CHROMS[1], '.psam'))

if(nrow(psam) > 100){
  target_qc_snplist<-plink_qc_snplist(pfile = opt$target_plink_chr_subset, chr = CHROMS, plink2 = opt$plink2, geno = opt$geno, maf = opt$maf, hwe = opt$hwe)
} else {
  target_qc_snplist<-plink_qc_snplist(pfile = opt$target_plink_chr_subset, chr = CHROMS, plink2 = opt$plink2, geno = opt$geno)
  log_add(log_file = log_file, message = 'Target sample size is <100 so only checking genotype missingness.')
}

###########
# QC reference
###########

ref_qc_snplist<-plink_qc_snplist(pfile = opt$ref_plink_chr_subset, chr = CHROMS, plink2 = opt$plink2, geno = opt$geno, maf = opt$maf, hwe = opt$hwe)

###########
# Harmonise target and reference genetic data
###########

# read in target pvar file
targ_pvar<-read_pvar(opt$target_plink_chr_subset, chr = CHROMS)

# read in reference pvar file
ref_pvar<-read_pvar(opt$ref_plink_chr_subset, chr = CHROMS)

# retain variants surviving QC
targ_pvar<-targ_pvar[targ_pvar$SNP %in% intersect(target_qc_snplist, ref_qc_snplist), ]
ref_pvar<-ref_pvar[ref_pvar$SNP %in% intersect(target_qc_snplist, ref_qc_snplist), ]

# insert IUPAC codes
targ_pvar$IUPAC <- snp_iupac(targ_pvar$A1, targ_pvar$A2)
ref_pvar$IUPAC <- snp_iupac(ref_pvar$A1, ref_pvar$A2)

# Identify SNPs present in both samples (allowing for strand flips)
# Identify SNPs that need to be flipped
target_ref<-merge(targ_pvar, ref_pvar, by='SNP')
flip <- detect_strand_flip(target_ref$IUPAC.x, target_ref$IUPAC.y)

flip_snplist<-NULL
if(sum(flip) > 0){
  flip_snplist<-target_ref$SNP.y[flip]
  log_add(log_file = log_file, message = paste0(sum(flip), 'variants will be flipped.'))
}

# Remove variants where IUPAC codes do not match (allowing for strand flips)
matched <- which((target_ref$IUPAC.x == target_ref$IUPAC.y) | flip)
target_ref<-target_ref[matched,]

log_add(log_file = log_file, message = paste0(nrow(target_ref),' variants match between target and reference after QC.'))

###########
# Identify list of LD independent SNPs
###########

log_add(log_file = log_file, message = 'Identifying LD independent SNPs based on reference data.')

# Subset ref_pvar to contain QC'd variants
ref_pvar<-ref_pvar[ref_pvar$SNP %in% target_ref$SNP,]

# Remove regions of high LD
ref_pvar <- remove_regions(dat = ref_pvar, regions = long_ld_coord)
log_add(log_file = log_file, message = paste0(nrow(ref_pvar),' variants after removal of LD high regions.'))

# Perform LD pruning
ld_indep <- plink_prune(pfile = opt$ref_plink_chr_subset, plink2 = opt$plink2, extract = ref_pvar$SNP, chr = CHROMS)
log_add(log_file = log_file, message = paste0(length(ld_indep),' independent variants retained.'))

###########
# Perform PCA based on reference
###########

log_add(log_file = log_file, message = 'Performing PCA based on reference.')

snp_weights<-plink_pca(pfile = opt$ref_plink_chr_subset, chr = CHROMS, plink2 = opt$plink2, extract = ld_indep, flip = flip_snplist, n_pc = opt$n_pcs)
fwrite(snp_weights, paste0(tmp_dir,'/ref.eigenvec.var'), row.names = F, quote=F, sep=' ', na='NA')

###
# Calculate PCs in the reference sample for scaling the target sample factor scores.
###

log_add(log_file = log_file, message = 'Computing reference PCs.')

# Calculate PCs in the reference
ref_pcs<-plink_score(pfile = opt$ref_plink_chr_subset, chr = CHROMS, plink2 = opt$plink2, score = paste0(tmp_dir,'/ref.eigenvec.var'), center = F)

# Scale across all individuals
ref_pcs_centre_scale <- score_mean_sd(scores = ref_pcs)

###
# Create model predicting reference populations
###

log_add(log_file = log_file, message = 'Deriving model predicting ref_pop groups.')

# Read in reference pop_data
pop_data <- read_pop_data(opt$pop_data)

# Scale the reference PCs
ref_pcs_scaled<-score_scale(score = ref_pcs, ref_scale = ref_pcs_centre_scale)
ref_pcs_scaled_pop<-merge(ref_pcs_scaled, pop_data, by=c('FID','IID'))

# Build model
model <- train(y=as.factor(ref_pcs_scaled_pop$POP), x=ref_pcs_scaled_pop[, grepl('PC',names(ref_pcs_scaled_pop)), with=F], method=opt$model_method, metric='logLoss', trControl=trainControl(method="cv", number=5, classProbs= TRUE, savePredictions = 'final', summaryFunction = multiClassSummary))

saveRDS(model$finalModel, paste0(opt$output,'.model.rds'))

#####
# Calculate PCs in target sample
#####

log_add(log_file = log_file, message = 'Calculating PCs in the target sample.')
targ_pcs<-plink_score(pfile = opt$target_plink_chr_subset, chr = CHROMS, plink2 = opt$plink2, score = paste0(tmp_dir,'/ref.eigenvec.var'), center = F)
targ_pcs_scaled<-score_scale(score = targ_pcs, ref_scale = ref_pcs_centre_scale)

###
# Create plot PC scores of target sample compared to the reference
###

log_add(log_file = log_file, message = 'Plotting target sample PCs on reference.')

# Combine ref and targ PCs
targ_pcs_scaled$POP<-'Target'
ref_pcs_targ_pcs<-rbind(ref_pcs_scaled_pop,targ_pcs_scaled)

PC_1_2<-ggplot(ref_pcs_targ_pcs[ref_pcs_targ_pcs$POP != 'Target',], aes(x=PC1,y=PC2, colour=POP)) +
  geom_point() +
	geom_point(data=ref_pcs_targ_pcs[ref_pcs_targ_pcs$POP == 'Target',], aes(x=PC1,y=PC2), colour='black', shape=21) +
	labs(title = "PCs 1 and 2", colour="") +
  theme_half_open() +
  background_grid() +
  theme(plot.title = element_text(hjust = 0.5))

PC_3_4<-ggplot(ref_pcs_targ_pcs[ref_pcs_targ_pcs$POP != 'Target',], aes(x=PC3,y=PC4, colour=POP)) +
  geom_point() +
	geom_point(data=ref_pcs_targ_pcs[ref_pcs_targ_pcs$POP == 'Target',], aes(x=PC3,y=PC4), colour='black', shape=21) +
	labs(title = "PCs 3 and 4", colour="") +
  theme_half_open() +
  background_grid() +
  theme(plot.title = element_text(hjust = 0.5))

png(paste0(opt$output,'.pc_plot.png'), units='px', res=300, width=4000, height=2000)
  plot_grid(PC_1_2,PC_3_4)
dev.off()

###
# Estimate probability of outcomes in model
###

log_add(log_file = log_file, message = 'Inferring population membership in target.')

# Read in model
pop_model_pred<-predict(object = model$finalModel, newx = data.matrix(targ_pcs_scaled[, grepl('PC',names(targ_pcs_scaled)), with=F]), type = "response", s=model$finalModel$lambdaOpt)
pop_model_pred<-as.data.frame.table(pop_model_pred)
pop_model_pred<-data.table(	FID=targ_pcs_scaled$FID,
														IID=targ_pcs_scaled$IID,
														pop=as.character(pop_model_pred$Var2),
														prob=round(pop_model_pred$Freq,3))

pop_model_pred<-dcast.data.table(pop_model_pred, formula=FID + IID~pop, value.var = "prob")

fwrite(pop_model_pred, paste0(opt$output,'.model_pred'), sep='\t')

# Create keep files based on the results
dir.create(paste0(opt$out_dir,'/keep_files/model_based'), recursive = T)
if(!is.na(opt$prob_thresh)){
  pop_model_pred$max_prob<-apply(pop_model_pred[,-1:-2], 1, max)
  pop_model_pred<-pop_model_pred[pop_model_pred$max_prob > opt$prob_thresh,]
  pop_model_pred$max_prob<-NULL
}

N_group<-NULL
for(i in names(pop_model_pred[,-1:-2])){
	tmp_keep<-pop_model_pred[apply(pop_model_pred[,-1:-2], 1, function(x) x[i] == max(x)),1:2]
	N_group<-rbind(N_group, data.frame(Group=i, N=nrow(tmp_keep)))
	fwrite(tmp_keep, paste0(opt$out_dir,'/keep_files/model_based/',i,'.keep'), sep=' ', col.names=F)
}

N_group<-rbind(N_group, data.frame(Group='Unassigned', N=nrow(targ_pcs_scaled) - nrow(pop_model_pred)))

sink(file = log_file, append = T)
cat('----------\n')
cat('N per group based on model:\n')
print.data.frame(N_group, row.names = FALSE, quote = FALSE, right = FALSE)
cat('----------\n')
sink()

if(opt$sd_rule){
  dir.create(paste0(opt$out_dir,'/keep_files/sd_based'), recursive = T)
  N_group<-NULL
  for(pop_i in unique(pop_data$POP)){

    # Calculate scale of PCs within reference population
    ref_pcs_scaled_i <- score_mean_sd(scores = ref_pcs, keep = pop_data[pop_data$POP == pop_i, c('FID','IID'), with=F])

    # Scale the target PC based on reference mean and SD
    targ_pcs_scaled_i<-score_scale(score = targ_pcs, ref_scale = ref_pcs_scaled_i)

    # Identify individuals with PCs <3
    targ_pcs_scaled_i<-targ_pcs_scaled_i[!apply(targ_pcs_scaled_i[,-1:-2], 1, function(x) any(x > 3 | x < -3)),]

    N_group<-rbind(N_group, data.frame(Group=pop_i, N=nrow(targ_pcs_scaled_i)))

    # Save keep file of individuals that fit the population
    fwrite(targ_pcs_scaled_i[,1:2], paste0(opt$out_dir,'/keep_files/sd_based/',pop_i,'.keep'), col.names=F, sep='\t')
  }

  sink(file = log_file, append = T)
  cat('----------\n')
  cat('N per group based on SD rule:\n')
  print.data.frame(N_group, row.names = FALSE, quote = FALSE, right = FALSE)
  cat('----------\n')
  sink()
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
