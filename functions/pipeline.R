#!/usr/bin/Rscript

if (!require("data.table", quietly = TRUE)) {
  library(data.table)
}

# Read in PGS
read_pgs <- function(config, name = NULL, pgs_methods = NULL, gwas = NULL, pop = NULL){

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

  pgs <- list()
  for (name_i in target_list$name) {
    # Read in keep_list to determine populations available
    keep_list_i <- fread(paste0(outdir,'/',name_i,'/ancestry/keep_list.txt'))

    if(!is.null(pop)){
      if(any(!(pop %in% keep_list_i$POP))){
        stop(paste0('Requested pop are not present in ',name_i,' sample.'))
      }
      keep_list_i <- keep_list_i[keep_list_i$POP %in% pop,]
    }

    pgs[[name_i]] <- list()
    for (pop_i in keep_list_i$POP) {
      pgs[[name_i]][[pop_i]] <- list()
      for(score_i in 1:nrow(score_file_list)){
        gwas_i <- score_file_list$name[score_i]
        pgs_method_i <- score_file_list$method[score_i]
        if (is.null(pgs[[name_i]][[pop_i]][[gwas_i]])) {
          pgs[[name_i]][[pop_i]][[gwas_i]] <- list()
        }
        pgs[[name_i]][[pop_i]][[gwas_i]][[pgs_method_i]] <-
            fread(
              paste0(
                outdir, '/', name_i, '/pgs/', pop_i, '/', pgs_method_i, '/',  gwas_i, '/', name_i, '-', gwas_i, '-', pop_i, '.profiles'
              )
            )
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
  if(pgs_method %in% pgs_group_methods){
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

  # Use most stringent p-value threshold of 0.05 as pseudo
  if(pgs_method == 'ptclump'){
    pseudo_val <- '0_1'
  }

  # Pseudoval only methods
  if(pgs_method == 'sbayesr'){
    pseudo_val <- 'SBayesR'
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
      cat('Returning result for EUR target population.\n')
      target_pop <- 'EUR'
    } else if(!is.null(target_pop) && !(target_pop %in% gwas_list$population)){
      cat(paste0('target_pop ', target_pop,' is not present in gwas_group ', gwas, '.\n'))
      cat('Returning result for EUR target population.\n')
      target_pop <- 'EUR'
    }
    pseudo_val <- paste0('targ_', target_pop, '_weighted')
  }

  if(tlprs){
    if(!is.null(target_pop) && target_pop == 'TRANS'){
      cat('No pseudovalidation for TRANS target population available for xwing.\n')
      cat('Returning result for EUR target population.\n')
      target_pop <- 'EUR'
    } else if(!is.null(target_pop) && !(target_pop %in% gwas_list$population)){
      cat(paste0('target_pop ', target_pop,' is not present in gwas_group ', gwas, '.\n'))
      cat('Returning result for EUR target population.\n')
      target_pop <- 'EUR'
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

