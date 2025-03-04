#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="GWAS ID [required]"),
  make_option("--config", action="store", default=NA, type='character',
              help="config file [required]")
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
if(is.null(opt$gwas)){
  stop('--gwas must be specified.\n')
}

# Identify outdir from config file
outdir <- read_param(config = opt$config, 'outdir', return_obj = F)

# Identify resdir from config file
resdir <- read_param(config = opt$config, 'resdir', return_obj = F)

# Identify refdir from config file
refdir <- read_param(config = opt$config, 'refdir', return_obj = F)

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(outdir,"/reference/gwas_sumstat/",opt$gwas,'/magma/set_extractor.log')
log_header(log_file = log_file, opt = opt, script = 'set_extractor.R', start.time = start.time)

# Read in significant and independent gene sets
if(!file.exists(paste0(outdir,'/reference/gwas_sumstat/',opt$gwas,'/magma/sig_indep_sets.txt'))){
  log_add(log_file = log_file, message = 'No sets were FDR significant.')
} else {
  
  set_enrich<-read.table(paste0(outdir,'/reference/gwas_sumstat/',opt$gwas,'/magma/sig_indep_sets.txt'), header=F)$V1
  
  # Read in .gmt file
  gmt_file <- read_param(config = opt$config, 'gene_sets', return_obj = F)
  set_annot<-readLines(gmt_file)
  set_ids<-sapply(strsplit(set_annot, '\t'),"[[",1)
  
  # Subset .gmt to contain enriched sets
  set_annot<-set_annot[set_ids %in% set_enrich]
  
  # Save list of significant sets
  dir.create(paste0(outdir,'/reference/gwas_sumstat/',opt$gwas,'/magma/snplists'))
  
  # Read in MAGMA gene locations file
  annot<-readLines(paste0(resdir, '/data/magma/NCBI37.3.genes.annot'))[-1:-2]
  annot<-strsplit(annot, '\t')
  annot_ids<-sapply(annot,"[[",1)
  
  ref <- read_pvar(paste0(refdir, '/ref.chr'))
  ref_snps <- ref$SNP
  
  for(set_i in 1:length(set_annot)){
    set_annot_i<-unlist(strsplit(set_annot[[set_i]], '\t'))
    
    set_id<-set_annot_i[1]
    genes<-set_annot_i[-1:-2]
    
    annot_subset<-annot[annot_ids %in% genes]
    
    snps<-unique(do.call(c, annot_subset))
    snps<-snps[grepl('^rs', snps)]
    snps<-snps[!is.na(snps)]
    
    snps<-snps[snps %in% ref_snps]
    
    write.table(snps, paste0(outdir,'/reference/gwas_sumstat/',opt$gwas,'/magma/snplists/',set_id,'.snplist'), col.names=F, row.names=F, quote=F)
  }
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()