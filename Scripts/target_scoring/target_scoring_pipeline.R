#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
library("optparse")

option_list = list(
make_option("--config", action="store", default=NULL, type='character',
    help="Pipeline configuration file [required]"),
make_option("--name", action="store", default=NULL, type='character',
    help="Name of target sample [required]"),
make_option("--population", action="store", default=NULL, type='character',
    help="Population in target sample to extract [required]"),
make_option("--plink2", action="store", default='plink2', type='character',
		help="Path PLINK v2 software binary [optional]"),
make_option("--gwas", action="store", default=NULL, type='character',
    help="Comma seperated list of names within gwas_list/score_list [required]"),
make_option("--pgs_method", action="store", default=NULL, type='character',
    help="Comma seperated list of polygenic scoring methods [required]"),
make_option("--n_cores", action="store", default=1, type='numeric',
    help="Number of cores to use [optional]"),
make_option("--test", action="store", default=NA, type='character',
    help="Specify number of SNPs to include [optional]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="Memory limit [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
library(GenoUtils)
library(data.table)
source('../functions/misc.R')
source_all('../functions')

# Check required inputs
if(is.null(opt$config)){
  stop('--config must be specified.\n')
}
if(is.null(opt$name)){
  stop('--name must be specified.\n')
}
if(is.null(opt$population)){
  stop('--population must be specified.\n')
}
if(is.null(opt$gwas)){
  stop('--gwas must be specified.\n')
}
if(is.null(opt$pgs_method)){
  stop('--pgs_method must be specified.\n')
}

# Read in outdir
outdir <- read_param(config = opt$config, param = 'outdir', return_obj = F)

# Create output directory
opt$output_dir <- paste0(outdir, '/', opt$name, '/pgs/', opt$population)
system(paste0('mkdir -p ',opt$output_dir))

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'target_scoring_pipeline.R', start.time = start.time)

# If testing, change CHROMS to chr value
if(!is.na(opt$test) && opt$test == 'NA'){
  opt$test<-NA
}
if(!is.na(opt$test)){
  CHROMS <- as.numeric(gsub('chr','',opt$test))
}

# Split opt$gwas
opt$gwas<-unlist(strsplit(gsub(' ','', opt$gwas), ','))
# Split opt$pgs_method
opt$pgs_method<-unlist(strsplit(gsub(' ','', opt$pgs_method), ','))

#####
# Combine score files
#####
# They are already sorted, have the same number of variants and have same A1
# Just append without duplicating the SNP, A1, and A2 columns, and insert method and gwas id in column names
first <- T
scale_files<-list()
for(pgs_method_i in opt$pgs_method){
  for(gwas_i in opt$gwas){
    if(first){
      system(paste0('zcat ', outdir, '/reference/pgs_score_files/', pgs_method_i,'/', gwas_i,'/ref-',gwas_i,".score.gz | sed '1 s/SCORE/",paste0(pgs_method_i,'_',gwas_i),"/' > ", tmp_dir,'/all_score.txt'))
      first <- F
    } else {
      system(paste0('zcat ', outdir, '/reference/pgs_score_files/', pgs_method_i,'/', gwas_i,'/ref-',gwas_i,".score.gz | cut -d' ' --complement -f1-3 - | sed '1 s/SCORE/",paste0(pgs_method_i,'_',gwas_i),"/' > ", tmp_dir,'/tmp_score.txt'))
      system(paste0("paste -d ' ' ", tmp_dir, '/all_score.txt ', tmp_dir, '/tmp_score.txt > ',tmp_dir, '/temp.txt && mv ', tmp_dir, '/temp.txt ', tmp_dir, '/all_score.txt'))
    }
    # Read in scale file and update Param
    scale_files[[paste0(pgs_method_i,'-',gwas_i)]]<-fread(paste0(outdir, '/reference/pgs_score_files/', pgs_method_i,'/', gwas_i,'/ref-',gwas_i,'-', opt$population,'.scale'))
    scale_files[[paste0(pgs_method_i,'-',gwas_i)]]$Param<-gsub('SCORE', paste0(pgs_method_i,'.',gwas_i), scale_files[[paste0(pgs_method_i,'-',gwas_i)]]$Param)
  }
}

# They will all be on the same strand and have the same variant IDs
# We just need to match the effect allele

# Create ref for all score to be matched to
refdir <- read_param(config = opt$config, param = 'refdir', return_obj = F)
ref<-read_pvar(paste0(refdir, '/ref.chr'), chr=CHROMS)
ref<-ref[, c('SNP','A1','A2'), with=F]

score_files<-list()
scale_files<-list()
for(pgs_method_i in opt$pgs_method){
  for(gwas_i in opt$gwas){
    # Read in score file
    tmp<-fread(paste0(outdir, '/reference/pgs_score_files/', pgs_method_i,'/', gwas_i,'/ref-',gwas_i,'.score.gz'))

    # Flip effects to match reference alleles
    tmp <- merge(ref, tmp, by = 'SNP', all.x=T, sort = F)
    flip <- which(tmp$A1.x != tmp$A1.y)
    tmp <- as.matrix(tmp[, -1:-5, drop = FALSE])
    tmp[flip, ] <- -tmp[flip, ,drop=F]
    score_files[[paste0(pgs_method_i,'-',gwas_i)]] <- tmp

    # Read in scale file
    scale_files[[paste0(pgs_method_i,'-',gwas_i)]]<-fread(paste0(outdir, '/reference/pgs_score_files/', pgs_method_i,'/', gwas_i,'/ref-',gwas_i,'-', opt$population,'.scale'))

    # Update names of scores to avoid conflicts between PGS methods
    # Could handle this when creating the score files
    colnames(score_files[[paste0(pgs_method_i,'-',gwas_i)]])<-gsub('SCORE', paste0(pgs_method_i,'.',gwas_i), colnames(score_files[[paste0(pgs_method_i,'-',gwas_i)]]))
    scale_files[[paste0(pgs_method_i,'-',gwas_i)]]$Param<-gsub('SCORE', paste0(pgs_method_i,'.',gwas_i), scale_files[[paste0(pgs_method_i,'-',gwas_i)]]$Param)
  }
}

# Concatenate score files
all_score <- cbind(ref, do.call(cbind, score_files))
rm(score_files)
rm(ref)

# Set all NA values to 0
all_score[is.na(all_score)]<-0

# Remove variants that have no effect
all_score <- all_score[apply(all_score[,-1:-3], 1, function(x) any(as.numeric(x) != 0)), ]

fwrite(all_score, paste0(tmp_dir,'/all_score.txt'), sep=' ', quote=F)

# Concatenate scale files
all_scale<-do.call(rbind, scale_files)

#####
# Perform polygenic risk scoring
#####

# Read in target_list
target_list <- read_param(config = opt$config, param = 'target_list', return_obj = T)

# Set params for plink_score
opt$target_plink_chr <- paste0(outdir, '/', opt$name, '/geno/', opt$name, '.ref.chr')
opt$target_keep <- paste0(outdir, '/', opt$name, '/ancestry/keep_files/model_based/', opt$population, '.keep')
opt$ref_freq_chr <- paste0(refdir, '/freq_files/', opt$population,'/ref.', opt$population,'.chr')

# Create a file listing SNPs to extract
fwrite(all_score[, 'SNP', with=F], paste0(tmp_dir,'/extract.txt'), col.names=F, quote=F)

log_add(log_file = log_file, message = 'Calculating polygenic scores in the target sample.')
scores <-
  plink_score(
    pfile = opt$target_plink_chr,
    chr = CHROMS,
    plink2 = opt$plink2,
    score = paste0(tmp_dir,'/all_score.txt'),
    extract = all_score$SNP,
    keep = opt$target_keep,
    frq = opt$ref_freq_chr,
    threads = opt$n_cores
  )

###
# Scale the polygenic scores based on the reference
###

log_add(log_file = log_file, message = 'Scaling target polygenic scores to the reference.')
scores_scaled<-score_scale(score=scores, ref_scale=all_scale)

###
# Write out the target sample scores
###

for(pgs_method_i in opt$pgs_method){
  for(gwas_i in opt$gwas){
    scores_scaled_i <- scores_scaled[, c('FID','IID', names(scores_scaled)[grepl(paste0(pgs_method_i,'.',gwas_i), names(scores_scaled))]), with=F]
    names(scores_scaled_i) <- gsub(paste0('^', pgs_method_i,'\\.'),'', names(scores_scaled_i))
    fwrite(scores_scaled_i, paste0(outdir, '/', opt$name,'/pgs/', opt$population,'/', pgs_method_i,'/', gwas_i,'/', opt$name,'-', gwas_i,'-',opt$population,'.profiles'), sep=' ', na='NA', quote=F)
  }
}

log_add(log_file = log_file, message = paste0('Saved polygenic scores.'))

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$output,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
