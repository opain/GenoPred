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

# Create temp directory
tmp_dir<-tempdir()

# Initiate log file
log_file <- paste0(outdir,"/reference/gwas_sumstat/",opt$gwas,'/magma/magma_set_conditional.log')
log_header(log_file = log_file, opt = opt, script = 'magma_set_conditional.R', start.time = start.time)

# Read in the MAGMA gene set enrichment results
sets_enrich<-fread(cmd=paste0("grep -v '^#' ",outdir,"/reference/gwas_sumstat/",opt$gwas,'/magma/magma_set_level.gsa.out'))

# Insert FULL_NAME column if not present
if(all(names(sets_enrich) != 'FULL_NAME')){
  sets_enrich$FULL_NAME<-sets_enrich$VARIABLE
}

# Remove gene sets with <5 genes present
sets_enrich<-sets_enrich[sets_enrich$NGENES >= 5,]

# Select FDR significant sets
# Note. this could be a parameter to be tuned
sets_enrich$P.FDR<-p.adjust(sets_enrich$P, method = 'fdr')
sets_enrich<-sets_enrich[sets_enrich$P.FDR <= 0.05,]

log_add(log_file = log_file, message = paste0(nrow(sets_enrich)," sets are FDR significant."))

# If more than 1 sig set, perform conditional analysis
if(nrow(sets_enrich) > 1){
  
  log_add(log_file = log_file, message = "Performing conditional analysis...")
  
  # Read in .gmt file
  gmt_file <- read_param(config = opt$config, 'gene_sets', return_obj = F)
  set_annot<-readLines(gmt_file)
  set_ids<-sapply(strsplit(set_annot, '\t'),"[[",1)
  
  # Subset .gmt to contain enriched sets
  set_annot<-set_annot[set_ids %in% sets_enrich$FULL_NAME]
  writeLines(set_annot, paste0(tmp_dir, "/sig_sets.gmt"))
  
  # Sort results by p-value
  sets_enrich<-sets_enrich[order(sets_enrich$P),]
  
  # Now condition each set on the most significant sets until all are independently significant
  i<-1
  
  sets_indep<-sets_enrich
  while(1){
    if(nrow(sets_indep) <= i){
      break
    }
    
    set_i<-sets_indep$FULL_NAME[1:i]
    
    log<-system(paste0(
      resdir, "/software/magma/magma",
      " --gene-results ",outdir,"/reference/gwas_sumstat/",opt$gwas,"/magma/magma_gene_level.genes.raw",
      " --set-annot ",tmp_dir, "/sig_sets.gmt",
      " --model direction-sets=greater condition-hide=",paste(set_i,collapse=','),
      " --out ",tmp_dir, "/res"
    ), intern = T)
    
    if(any(grepl('ERROR - running gene-level regression: could not invert design matrix of conditioned-on variables; variables are collinear with each other', log))){
      print(log)
      print('ERROR: There was too much multicolinearity between sets.')
      q()
    }
    
    # Read in the results
    cond_res<-fread(cmd=paste0("grep -v '^#' ",tmp_dir, '/res.gsa.out'))
    
    # Insert FULL_NAME column if not present
    if(all(names(cond_res) != 'FULL_NAME')){
      cond_res$FULL_NAME<-cond_res$VARIABLE
    }
    
    # Remove sets from sets_indep that are no longer significant (P>0.01)
    cond_res<-cond_res[cond_res$P < 0.01,]
    sets_indep<-sets_indep[sets_indep$FULL_NAME %in% c(set_i,cond_res$FULL_NAME),]
    
    i<-i+1
  }
  
  # Save file listing significant and independent sets
  write.table(sets_indep$FULL_NAME, paste0(outdir,"/reference/gwas_sumstat/",opt$gwas,'/magma/sig_indep_sets.txt'), row.names=F, col.names=F, quote=F)
  
  log_add(log_file = log_file, message = paste0(nrow(sets_indep), " independent sets remain."))
}

# If 1 sig set, no conditional analysis required
if(nrow(sets_enrich) == 1){
  log_add(log_file = log_file, message = "No conditional analysis required.")
  
  # Save file listing significant and independent sets
  write.table(sets_enrich$FULL_NAME, paste0(outdir,"/reference/gwas_sumstat/",opt$gwas,'/magma/sig_indep_sets.txt'), row.names=F, col.names=F, quote=F)
}

if(nrow(sets_enrich) == 0){
  log_add(log_file = log_file, message = 'No conditional analysis required')
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = log_file, append = T)
cat('Analysis finished at', as.character(end.time),'\n')
cat('Analysis duration was', as.character(round(time.taken,2)), attr(time.taken, 'units'), '\n')
sink()