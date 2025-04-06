#!/usr/bin/Rscript

if (!require("data.table", quietly = TRUE)) {
  library(data.table)
}

# Read in PGS
read_pgs <- function(config, name = NULL, pgs_methods = NULL, gwas = NULL, pop = NULL, pseudo_only = F){
  # Read in target_list
  target_list <- read_param(config = config, param = 'target_list')
  if(!is.null(name)){
    if(any(!(name %in% target_list$name))){
      stop('Requested target samples are not present in target_list')
    }
    name_i <- name
    target_list <- target_list[target_list$name %in% name_i,]
  }

  # Identify score files
  score_file_list <- list_score_files(config)

  # Subset requested gwas
  if(!is.null(gwas)){
    if(any(!(gwas %in% score_file_list$name))){
      stop('Requested GWAS are not present in gwas_list/score_list')
    }
    score_file_list<-score_file_list[score_file_list$name %in% gwas,]
  }

  # Subset requested pgs_methods
  if(!is.null(pgs_methods)){
    if(any(!(pgs_methods %in% score_file_list$method))){
      stop('Requested PGS methods are not present in gwas_list/score_list')
    }
    score_file_list<-score_file_list[score_file_list$method %in% pgs_methods,]
  }

  # Identify outdir parameter
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)
  
  # Identify pgs_scaling parameter
  pgs_scaling <- read_param(config = config, param = 'pgs_scaling', return_obj = F)
  
  pgs <- list()
  for (name_i in target_list$name) {
    pops<-NULL
    if('continuous' %in% pgs_scaling){
      pops <- c('TRANS', pops)
    }
    if('discrete' %in% pgs_scaling){
      # Read in keep_list to determine populations available
      keep_list_i <- fread(paste0(outdir,'/',name_i,'/ancestry/keep_list.txt'))
      pops <- c(pops, keep_list_i$POP)
    }
    if(!is.null(pop)){
      if(!any('discrete' %in% pgs_scaling) & any(pop != 'TRANS')){
        stop(paste0('Requested pop are not present in ',name_i,' sample. Only PGS adjusted using continuous ancestry correction are available due to pgs_scaling parameter in configfile.'))
      }
      if(any(!(pop %in% pops))){
        stop(paste0('Requested pop are not present in ',name_i,' sample.'))
      }
      pops <- pops[pops %in% pop]
    }

    pgs[[name_i]] <- list()
    for (pop_i in pops) {
      pgs[[name_i]][[pop_i]] <- list()
      for(score_i in 1:nrow(score_file_list)){
        gwas_i <- score_file_list$name[score_i]
        pgs_method_i <- score_file_list$method[score_i]
        if (is.null(pgs[[name_i]][[pop_i]][[gwas_i]])) {
          pgs[[name_i]][[pop_i]][[gwas_i]] <- list()
        }
        file_i<-paste0(outdir, '/', name_i, '/pgs/', pop_i, '/', pgs_method_i, '/',  gwas_i, '/', name_i, '-', gwas_i, '-', pop_i, '.profiles')
        if(pseudo_only){
          pseudo_param <- find_pseudo(config = config, gwas = gwas_i, target_pop = pop_i, pgs_method = pgs_method_i)

          score_header <-
            fread(file_i, nrows = 1)
          score_cols <-
            which(names(score_header) %in% c('FID', 'IID', paste0(gwas_i, '_',pseudo_param)))
          
          pgs[[name_i]][[pop_i]][[gwas_i]][[pgs_method_i]] <-
            fread(cmd = paste0("cut -d' ' -f ", paste0(score_cols, collapse=','), " ", file_i))
        } else {
          pgs[[name_i]][[pop_i]][[gwas_i]][[pgs_method_i]] <- fread(file_i)
        }
      }
    }
  }

  return(pgs)
}

# Read in PGS
read_pgs_2 <- function(config, name = NULL, pgs_methods = NULL, gwas = NULL, pop = NULL, pseudo_only = F, partitioned = F){
  # Identify outdir parameter
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)
  
  # Read in target_list
  target_list <- read_param(config = config, param = 'target_list')
  if(!is.null(name)){
    if(any(!(name %in% target_list$name))){
      stop('Requested target samples are not present in target_list')
    }
    name_i <- name
    target_list <- target_list[target_list$name %in% name_i,]
  }
  
  # Identify score files
  score_file_list <- list_score_files(config)
  
  # If partitioned, restrict to single source methods, and gwas with sig sets
  if(partitioned){
    score_file_list <- score_file_list[!(score_file_list$method %in% pgs_group_methods) & !grepl('tlprs|leopard', score_file_list$method),]
    
    set_reporter_file <- paste0(outdir, '/reference/gwas_sumstat/set_reporter.txt')
    set_reporter<-fread(set_reporter_file)
    score_file_list<-score_file_list[score_file_list$name %in% set_reporter$name[set_reporter$n_sig > 0],]
    
    part<-'.partitioned'
  }
  
  # Subset requested gwas
  if(!is.null(gwas)){
    if(any(!(gwas %in% score_file_list$name))){
      stop('Requested GWAS are not present in gwas_list/score_list')
    }
    score_file_list<-score_file_list[score_file_list$name %in% gwas,]
  }
  
  # Subset requested pgs_methods
  if(!is.null(pgs_methods)){
    if(any(!(pgs_methods %in% score_file_list$method))){
      stop('Requested PGS methods are not present in gwas_list/score_list')
    }
    score_file_list<-score_file_list[score_file_list$method %in% pgs_methods,]
  }
  
  # Identify pgs_scaling parameter
  pgs_scaling <- read_param(config = config, param = 'pgs_scaling', return_obj = F)
  
  pgs <- list()
  for (name_i in target_list$name) {
    pops<-NULL
    if('continuous' %in% pgs_scaling){
      pops <- c('TRANS', pops)
    }
    if('discrete' %in% pgs_scaling){
      # Read in keep_list to determine populations available
      keep_list_i <- fread(paste0(outdir,'/',name_i,'/ancestry/keep_list.txt'))
      pops <- c(pops, keep_list_i$POP)
    }
    if(!is.null(pop)){
      if(!any('discrete' %in% pgs_scaling) & any(pop != 'TRANS')){
        stop(paste0('Requested pop are not present in ',name_i,' sample. Only PGS adjusted using continuous ancestry correction are available due to pgs_scaling parameter in configfile.'))
      }
      if(any(!(pop %in% pops))){
        stop(paste0('Requested pop are not present in ',name_i,' sample.'))
      }
      pops <- pops[pops %in% pop]
    }
    
    pgs[[name_i]] <- list()
    for (pop_i in pops) {
      pgs[[name_i]][[pop_i]] <- list()
      for(score_i in 1:nrow(score_file_list)){
        gwas_i <- score_file_list$name[score_i]
        pgs_method_i <- score_file_list$method[score_i]
        if (is.null(pgs[[name_i]][[pop_i]][[gwas_i]])) {
          pgs[[name_i]][[pop_i]][[gwas_i]] <- list()
        }
        file_i<-paste0(outdir, '/', name_i, '/pgs/', pop_i, '/', pgs_method_i, '/',  gwas_i, '/', name_i, '-', gwas_i, '-', pop_i, part, '.profiles')
        if(pseudo_only){
          pseudo_param <- find_pseudo(config = config, gwas = gwas_i, target_pop = pop_i, pgs_method = pgs_method_i)
          
          score_header <-
            fread(file_i, nrows = 1)
          score_cols <-
            which(names(score_header) %in% c('FID', 'IID', paste0(gwas_i, '_',pseudo_param)))
          
          pgs[[name_i]][[pop_i]][[gwas_i]][[pgs_method_i]] <-
            fread(cmd = paste0("cut -d' ' -f ", paste0(score_cols, collapse=','), " ", file_i))
        } else {
          pgs[[name_i]][[pop_i]][[gwas_i]][[pgs_method_i]] <- fread(file_i)
        }
      }
    }
  }
  
  return(pgs)
}

# Create function to read in parameters in the config file
read_param <- function(config, param, return_obj = T){
  library(yaml)

  # Read in the config file
  config_file <- read_yaml(config)

  if(all(names(config_file) != param)){
    # Check default config file
    config_file <- read_yaml('config.yaml')

    if(all(names(config_file) != param)){
      cat(param, 'parameter is not present in user specified config file or default config file.\n')
      return(NULL)
    } else {
      cat(param, 'parameter is not present in user specified config file, so will use value in default config file.\n')
    }
  }

  # Identify value for param
  file <- config_file[[param]]
  file[file == 'NA']<-NA

  # If resdir, and NA, set to 'resources'
  if(param == 'resdir'){
    if(is.na(file)){
      file <- 'resources'
    }
  }

  # If refdir, and NA, set to '<resdir>/data/ref'
  if(param == 'refdir'){
    if(is.na(file)){
      resdir <- read_param(config = config, param = 'resdir', return_obj = F)
      file <- paste0(resdir, '/data/ref')
    }
  }

  if(return_obj){
    if(!is.na(file)){
      obj <- fread(file)
    } else {
      obj <- NULL
    }
    return(obj)
  } else {
    file <- file[order(file)]
    return(file)
  }
}

# Read in ancestry classifications
read_ancestry <- function(config, name){

  # Read in the config file
  config_file <- readLines(config)

  # Identify outdir
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  # Read keep_list
  keep_list <- fread(paste0(outdir,'/',name,'/ancestry/keep_list.txt'))

  # Read in keep lists
  keep_files <- list()
  for(pop in keep_list$POP){
    keep_files[[pop]] <- fread(keep_list$file[keep_list$POP == pop], header = F)
  }

  # Read in model predictions
  model_pred <- fread(paste0(outdir,'/',name,'/ancestry/',name,'.Ancestry.model_pred'))

  # Read in ancestry log file
  log <- readLines(paste0(outdir,'/',name,'/ancestry/',name,'.Ancestry.log'))

  output <- list(
    keep_list = keep_list,
    keep_files = keep_files,
    model_pred = model_pred,
    log = log
  )

  return(output)
}

# Return score corresponding to pseudovalidation
find_pseudo <- function(config, gwas, pgs_method, target_pop = NULL){

  if(length(pgs_method) > 1){
    stop('Only one pgs_method can be specified at a time')
  }
  if(length(gwas) > 1){
    stop('Only one gwas can be specified at a time')
  }
  if(length(target_pop) > 1){
    stop('Only one target_pop can be specified at a time')
  }
  if(pgs_method %in% pgs_group_methods & is.null(target_pop)){
    stop('target_pop must be specified when using multi-ancestry PGS method')
  }

  # Read in gwas_list
  gwas_list <- read_param(config = config, param = 'gwas_list')

  # Read in gwas_groups
  gwas_groups <- read_param(config = config, param = 'gwas_groups')

  # If pgs_method is multi-source, subset gwas_list to gwas in relevant group
  if(grepl(paste0('^', pgs_group_methods, collapse = '|'), pgs_method)){
    gwas_list <- gwas_list[gwas_list$name %in% unlist(strsplit(gwas_groups$gwas[gwas_groups$name == gwas], ','))]
  }

  # Identify score files
  score_file_list <- list_score_files(config)

  # Subset requested gwas
  if(!is.null(gwas)){
    if(any(!(gwas %in% score_file_list$name))){
      stop('Requested GWAS are not present in gwas_list/score_list')
    }
    score_file_list<-score_file_list[score_file_list$name %in% gwas,]
  }

  # Subset requested pgs_methods
  if(!is.null(pgs_method)){
    if(any(!(pgs_method %in% score_file_list$method))){
      stop('Requested PGS method are not present in gwas_list/score_list')
    }
    score_file_list<-score_file_list[score_file_list$method %in% pgs_method,]
  }

  # Find outdir
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  # If TLPRS, find pseudo param, and then edit value for TLPRS
  tlprs <- ifelse(grepl('tlprs', pgs_method), T, F)
  pgs_method <- gsub('tlprs_', '', pgs_method)
  if(tlprs && pgs_method %in% c('lassosum','megaprs')){
    if(!is.null(target_pop) && target_pop == 'TRANS'){
      cat('No pseudovalidation for TRANS target population available for ', pgs_method, '\n')
      cat(paste0('Returning result for ', gwas_list$population[1],' target population.\n'))
      target_pop <- gwas_list$population[1]
    }
    if(!is.null(target_pop) && target_pop %in% gwas_list$population){
      # Note. selecting pseudoval from non-target GWAS, as this the score file going into TLPRS
      gwas <- gwas_list$name[gwas_list$population != target_pop]
    } else {
      cat(paste0('target_pop ', target_pop,' is not present in gwas_group ', gwas, '.\n'))
      cat(paste0('Returning result for ', gwas_list$population[1],' target population.\n'))
      target_pop <- gwas_list$population[1]
      # Note. selecting pseudoval from non-target GWAS, as this the score file going into TLPRS
      gwas <- gwas_list$name[gwas_list$population != target_pop]
    }
  }

  # Use most stringent p-value threshold of 0.05 as pseudo
  if(pgs_method == 'ptclump'){
    pseudo_val <- '0_1'
  }

  # Pseudoval only methods
  if(pgs_method == 'sbayesr'){
    pseudo_val <- 'SBayesR'
  }
  if(pgs_method == 'sbayesrc'){
    pseudo_val <- 'SBayesRC'
  }
  if(pgs_method == 'quickprs'){
    pseudo_val <- 'quickprs'
  }

  # Retrieve pseudoval param
  if(pgs_method == 'dbslmm'){
    pseudo_val <- 'DBSLMM_1'
  }
  if(pgs_method == 'ldpred2'){
    pseudo_val <- 'beta_auto'
  }
  if(pgs_method == 'prscs'){
    pseudo_val <- 'phi_auto'
  }

  if(pgs_method == 'megaprs'){
    # Read in megaprs log file
    log <- readLines(paste0(outdir,'/reference/pgs_score_files/',pgs_method,'/',gwas,'/ref-',gwas,'.log'))
    log <- log[grepl('identified as the best with correlation', log)]
    pseudoval <- gsub(' .*','', gsub('Model ', '', log))
    pseudo_val <- paste0('ldak_Model', pseudoval)
  }
  if(pgs_method == 'lassosum'){
    # Read in megaprs log file
    log <- readLines(paste0(outdir,'/reference/pgs_score_files/',pgs_method,'/',gwas,'/ref-',gwas,'.log'))
    s_val <- gsub('.* ', '', log[grepl('^s = ', log)])
    lambda_val <- gsub('.* ', '', log[grepl('^lambda = ', log)])
    pseudo_val <- paste0('s', s_val, '_lambda', lambda_val)
  }

  # If pgs_method is external, return the only score
  if(pgs_method == 'external'){
    pseudo_val <- 'external'
  }

  # Multi-population methods
  if(pgs_method == 'prscsx'){
    pseudo_val <- 'META_phi_auto'
  }
  if(pgs_method == 'xwing'){
    if(!is.null(target_pop) && target_pop == 'TRANS'){
      cat('No pseudovalidation for TRANS target population available for xwing.\n')
      cat(paste0('Returning result for ', gwas_list$population[1],' target population.\n'))
      target_pop <- gwas_list$population[1]
    } else if(!is.null(target_pop) && !(target_pop %in% gwas_list$population)){
      cat(paste0('target_pop ', target_pop,' is not present in gwas_group ', gwas, '.\n'))
      cat(paste0('Returning result for ', gwas_list$population[1],' target population.\n'))
      target_pop <- gwas_list$population[1]
    }
    pseudo_val <- paste0('targ_', target_pop, '_weighted')
  }
  if(grepl('_multi$', pgs_method)){
    if(!is.null(target_pop) && target_pop == 'TRANS'){
      cat('No pseudovalidation for TRANS target population available for xwing.\n')
      cat(paste0('Returning result for ', gwas_list$population[1],' target population.\n'))
      target_pop <- gwas_list$population[1]
    } else if(!is.null(target_pop) && !(target_pop %in% gwas_list$population)){
      cat(paste0('target_pop ', target_pop,' is not present in gwas_group ', gwas, '.\n'))
      cat(paste0('Returning result for ', gwas_list$population[1],' target population.\n'))
      target_pop <- gwas_list$population[1]
    }
    pseudo_val <- paste0('targ_', target_pop, '_weighted')
  }

  if(tlprs){
    if(!is.null(target_pop) && target_pop == 'TRANS'){
      cat('No pseudovalidation for TRANS target population available for TLPRS\n')
      cat(paste0('Returning result for ', gwas_list$population[1],' target population.\n'))
      target_pop <- gwas_list$population[1]
    } else if(!is.null(target_pop) && !(target_pop %in% gwas_list$population)){
      cat(paste0('target_pop ', target_pop,' is not present in gwas_group ', gwas, '.\n'))
      cat(paste0('Returning result for ', gwas_list$population[1],' target population.\n'))
      target_pop <- gwas_list$population[1]
    }
    pseudo_val <- paste0('targ_', target_pop, '_', pseudo_val, '_TLPRS_61')
  }
  return(pseudo_val)
}

# Read in lassosum pseudoval results
read_pseudo_r <- function(config, gwas){

  if(length(gwas) > 1){
    stop('Only one gwas can be specified at a time')
  }

  # Find outdir param
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  # Read in lassosum log file
  log <- readLines(paste0(outdir,'/reference/pgs_score_files/lassosum/',gwas,'/ref-',gwas,'.log'))
  r <- as.numeric(gsub('value = ','',log[grepl('value = ', log)]))

  return(r)
}

# Read in reference PGS
# Read in TRANS scores (adjusted for ancestry), and restrict to pseudovalidated models
read_reference_pgs <- function(config){
  
  # Identify score files
  score_file_list <- list_score_files(config)
  
  # Identify outdir parameter
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)
  
  pgs <- list()
  for(score_i in 1:nrow(score_file_list)){
    gwas_i <- score_file_list$name[score_i]
    pgs_method_i <- score_file_list$method[score_i]
    if (is.null(pgs[[gwas_i]])) {
      pgs[[gwas_i]] <- list()
    }
    pgs[[gwas_i]][[pgs_method_i]] <-
      fread(
        paste0(
          outdir, '/reference/pgs_score_files/', pgs_method_i, '/',  gwas_i, '/ref-', gwas_i, '-TRANS.profiles'
        )
      )
    pseudo_param <- find_pseudo(config = config, gwas = gwas_i, pgs_method = pgs_method_i, target_pop = 'TRANS')
    pgs[[gwas_i]][[pgs_method_i]]<-pgs[[gwas_i]][[pgs_method_i]][,c('FID','IID',paste0('SCORE_',pseudo_param)), with=F]
  }
  
  return(pgs)
}

