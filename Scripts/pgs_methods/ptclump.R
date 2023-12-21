#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--ref_plink_chr", action="store", default=NULL, type='character',
		help="Path to per chromosome reference PLINK files [required]"),
make_option("--ref_keep", action="store", default=NULL, type='character',
    help="Keep file to subset individuals in reference for clumping [optional]"),
make_option("--pop_data", action="store", default=NULL, type='character',
		help="File containing the population code and location of the keep file [required]"),
make_option("--plink", action="store", default='plink', type='character',
    help="Path PLINKv1.9 software binary [required]"),
make_option("--plink2", action="store", default='plink2', type='character',
    help="Path PLINKv2 software binary [required]"),
make_option("--output", action="store", default='./Output', type='character',
		help="Path for output files [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]"),
make_option("--sumstats", action="store", default=NULL, type='character',
		help="GWAS summary statistics in LDSC format [required]"),
make_option("--pTs", action="store", default='1e-8,1e-6,1e-4,1e-2,0.1,0.2,0.3,0.4,0.5,1', type='character',
		help="List of p-value thresholds for scoring [optional]"),
make_option("--extract", action="store", default=NULL, type='character',
    help="File listing SNPs to extract for polygenic scoring [optional]"),
make_option("--nested", action="store", default=T, type='logical',
    help="Specify as F to use non-overlapping p-value intervals [optional]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--top_hla", action="store", default=T, type='logical',
		help="Retain only top assocaited variant in HLA/MHC region [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../Scripts/functions/misc.R')

# Check required inputs
if(is.null(opt$ref_plink_chr)){
  stop('--ref_plink_chr must be specified.\n')
}
if(is.null(opt$sumstats)){
  stop('--sumstats must be specified.\n')
}
if(is.null(opt$pop_data)){
  stop('--pop_data must be specified.\n')
}

# Create output directory
opt$output_dir<-paste0(dirname(opt$output),'/')
system(paste0('mkdir -p ',opt$output_dir))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'ptclump.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

#####
# Format pT option
#####

opt$pTs<-as.numeric(unlist(strsplit(opt$pTs,',')))

#####
# Read in sumstats
#####

log_add(log_file = log_file, message = 'Reading in GWAS.')  

# Read in, check and format GWAS summary statistics
gwas <- read_sumstats(sumstats = opt$sumstats, chr = CHROMS, log_file = log_file, extract = opt$extract, req_cols = c('CHR','BP','SNP','A1','A2','BETA','P'))

#####
# Prepare for ptclump
#####

if(opt$top_hla){
  # Assumes BP column is GRCh37
  hla <- gwas[(gwas$CHR == 6 & gwas$BP > 28e6 & gwas$BP < 34e6),]
  top_hla <- hla$SNP[hla$P == min(hla$P)][1]
  gwas <- gwas[!(gwas$CHR == 6 & gwas$BP > 28e6 & gwas$BP < 34e6 & gwas$SNP != top_hla),]
  log_add(log_file = log_file, message = c('Extracted top variant in HLA/MHC region.', paste0(nrow(gwas), ' variants remain.')))  
}

#####
# Clump SNPs in GWAS based on LD in the reference
#####

# Record start time for test
if(!is.na(opt$test)){
  test_start.time<-test_start(log_file = log_file)
}

# Peforming LD-based clumping
if(!is.null(opt$keep)){
  log_add(log_file = log_file, message = 'Restricting sample to keep file for clumping.')
}

clumped <- plink_clump(bfile = opt$ref_plink_chr, chr = CHROMS, sumstats = gwas, keep = opt$ref_keep, log_file = log_file)

#####
# Create score files
#####

log_add(log_file = log_file, message = 'Creating score file.')  

# Restrict GWAS to clumped variants
gwas <- gwas[(gwas$SNP %in% clumped),]

# Retain pTs with at least one variant
opt$pTs<-opt$pTs[opt$pTs > min(gwas$P)]

# Create range_list file based on specified p-value thresholds
if(opt$nested == T){
  range_list<-data.frame(	Name=paste0('S',1:length(opt$pTs)),
                          pT0=0,
                          pT1=opt$pTs)
} else {
  range_list<-data.frame(	Name=paste0('S',1:length(opt$pTs)),
                          pT0=c(0,opt$pTs[-length(opt$pTs)]),
                          pT1=opt$pTs)
}

score <- gwas[, c('SNP','A1','A2'), with=F]
for(pT in 1:nrow(range_list)){
  tmp <- gwas$BETA
  tmp[!(gwas$P > range_list$pT0[pT] & gwas$P < range_list$pT1[pT])] <- 0
  score[[paste0('SCORE_',range_list$pT0[pT],'_',range_list$pT1[pT])]] <- tmp
}

fwrite(score, paste0(opt$output,'.score'), col.names=T, sep=' ', quote=F)

if(file.exists(paste0(opt$output,'.score.gz'))){
  system(paste0('rm ',opt$output,'.score.gz'))
}
system(paste0('gzip ',opt$output,'.score'))

# Record end time of test
if(!is.na(opt$test)){
  test_finish(log_file = log_file, test_start.time = test_start.time)
}

###
# Count the number of SNPs to be included at each pT
###

for(i in 1:nrow(range_list)){
	range_list$NSNP[i]<-sum(gwas$P > range_list$pT0[i] & gwas$P < range_list$pT1[i])
}

fwrite(range_list, paste0(opt$output,'.NSNP_per_pT'), sep='\t')

####
# Calculate mean and sd of polygenic scores at each threshold
####

log_add(log_file = log_file, message = 'Calculating polygenic scores in reference.')  

# Calculate scores in the full reference
ref_pgs<-calc_score(bfile = opt$ref_plink_chr, chr = CHROMS, plink2 = opt$plink2, score = paste0(opt$output,'.score.gz'))

# Calculate scale within each reference population
pop_data<-fread(opt$pop_data)

for(pop_i in unique(pop_data$POP)){
  ref_pgs_scale_i <- score_mean_sd(scores = ref_pgs, keep = pop_data[pop_data$POP == pop_i, c('FID','IID'), with=F])
  fwrite(ref_pgs_scale_i, paste0(opt$output, '.', pop_i, '.scale'), row.names = F, quote=F, sep=' ', na='NA')
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
