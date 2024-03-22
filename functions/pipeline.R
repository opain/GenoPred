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

  # Read in gwas_list
  gwas_list <- read_param(config = config, param = 'gwas_list')

  # Read in score_list
  score_list <- read_param(config = config, param = 'score_list')

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  if(!is.null(score_list)){
    # Read in score_reporter output
    score_reporter <- fread(paste0(outdir, "/reference/pgs_score_files/external/score_report.txt"))
    score_list <- merge(score_list, score_reporter, by='name')

    # Remove scores that did not pass ref harmonisation
    score_list <- score_list[score_list$pass == T,]
  }

  if(!is.null(gwas)){
    if(!is.null(score_list)){
      full_gwas_list <- c(gwas_list$name, score_list$name)
    } else {
      full_gwas_list <- gwas_list$name
    }

    if(any(!(gwas %in% full_gwas_list))){
      stop('Requested GWAS are not present in gwas_list/score_list')
    }
    gwas_list <- gwas_list[gwas_list$name %in% gwas,]

    if(!is.null(score_list)){
      score_list <- score_list[score_list$name %in% gwas,]
    }
  }

  # Identify PGS methods to be included
  pgs_methods_list <- read_param(config = config, param = 'pgs_methods', return_obj = F)

  if(!is.null(pgs_methods)){
    if(!is.null(score_list)){
      if(any(!(pgs_methods %in% c(pgs_methods_list, 'external')))){
        stop('Requested pgs_methods are not present in pgs_methods in config')
      }
    } else {
      if(any(!(pgs_methods %in% pgs_methods_list))){
        stop('Requested pgs_methods are not present in pgs_methods in config')
      }
    }
    pgs_methods_list <- pgs_methods_list[pgs_methods_list %in% pgs_methods]
  }

  # Define PGS methods applied to non-EUR GWAS
  pgs_methods_noneur <- c('ptclump','lassosum','megaprs','prscs','dbslmm')

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
      for (gwas_i in gwas_list$name) {
        pgs[[name_i]][[pop_i]][[gwas_i]] <- list()

        for (pgs_method_i in pgs_methods_list) {
          if (gwas_list$population[gwas_list$name == gwas_i] == 'EUR' | (gwas_list$population[gwas_list$name == gwas_i] != 'EUR' & (pgs_method_i %in% pgs_methods_noneur))) {
            pgs[[name_i]][[pop_i]][[gwas_i]][[pgs_method_i]] <-
              fread(
                paste0(
                  outdir, '/', name_i, '/pgs/', pop_i, '/', pgs_method_i, '/',  gwas_i, '/', name_i, '-', gwas_i, '-', pop_i, '.profiles'
                )
              )
          }
        }
      }
      if(!is.null(score_list)){
        for (score_i in score_list$name) {
          pgs[[name_i]][[pop_i]][[score_i]] <- list()
          pgs_method_i <- 'external'
          pgs[[name_i]][[pop_i]][[score_i]][[pgs_method_i]] <-
            fread(
              paste0(
                outdir, '/', name_i, '/pgs/', pop_i, '/', pgs_method_i, '/',  score_i, '/', name_i, '-', score_i, '-', pop_i, '.profiles'
              )
            )
        }
      }
    }
  }

  return(pgs)
}

# Create function to read in parameters in the config file
read_param <- function(config, param, return_obj = T){

  # Read in the config file
  config_file <- readLines(config)

  if(all(grepl(paste0('^',param,':'), config_file) == F)){
    # Check default config file
    config_file <- readLines('config.yaml')

    if(all(grepl(paste0('^',param,':'), config_file) == F)){
      cat('Requested parameter is not present in user specified config file or default config file.')
      return(NULL)
    } else {
      cat('Parameter is not present in user specified config file, so will use value in default config file.')
    }
  }

  # Identify value for param
  file <- gsub(paste0(param,': '), '', config_file[grepl(paste0('^',param,':'), config_file)])
  file[file == 'NA']<-NA

  if(return_obj){
    if(!is.na(file)){
      obj <- fread(file)
    } else {
      obj <- NULL
    }
    return(obj)
  } else {

    file <- unlist(strsplit(gsub("'", '', gsub(']', '', gsub('\\[', '', file))),','))
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
find_pseudo <- function(config, gwas, pgs_method){

  if(length(pgs_method) > 1){
    stop('Only one pgs_method can be specified at a time')
  }
  if(length(gwas) > 1){
    stop('Only one gwas can be specified at a time')
  }

  # Read in gwas_list
  gwas_list <- read_param(config = config, param = 'gwas_list')

  # Read in score_list
  score_list <- read_param(config = config, param = 'score_list')

  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  if(!is.null(score_list)){
    # Read in score_reporter output
    score_reporter <- fread(paste0(outdir, "/reference/pgs_score_files/external/score_report.txt"))
    score_list <- merge(score_list, score_reporter, by='name')

    # Remove scores that did not pass ref harmonisation
    score_list <- score_list[score_list$pass == T,]
  }

  # Find outdir
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  if(!is.null(gwas)){
    if(!is.null(score_list)){
      full_gwas_list <- c(gwas_list$name, score_list$name)
    } else {
      full_gwas_list <- gwas_list$name
    }

    if(any(!(gwas %in% full_gwas_list))){
      stop('Requested GWAS are not present in gwas_list/score_list')
    }
    gwas_list <- gwas_list[gwas_list$name %in% gwas,]

    if(!is.null(score_list)){
      score_list <- score_list[score_list$name %in% gwas,]
    }
  }

  # Identify PGS methods to be included
  pgs_methods_list <- read_param(config = config, param = 'pgs_methods', return_obj = F)

  if(!is.null(pgs_method)){
    if(!is.null(score_list)){
      if(any(!(pgs_method %in% c(pgs_methods_list, 'external')))){
        stop('Requested pgs_method are not present in pgs_methods in config')
      }
    } else {
      if(any(!(pgs_method %in% pgs_methods_list))){
        stop('Requested pgs_method are not present in pgs_methods in config')
      }
    }
    pgs_methods_list <- pgs_methods_list[pgs_methods_list %in% pgs_method]
  }

  # Identify outdir parameter
  outdir <- read_param(config = config, param = 'outdir', return_obj = F)

  # Use most stringent p-value threshold of 0.05 as pseudo
  if(pgs_method == 'ptclump'){
    return('0_1')
  }

  # Pseudoval only methods
  if(pgs_method == 'sbayesr'){
    return('SBayesR')
  }


  # Retrieve pseudoval param
  if(pgs_method == 'dbslmm'){
    return('DBSLMM_1')
  }
  if(pgs_method == 'ldpred2'){
    return('beta_auto')
  }
  if(pgs_method == 'prscs'){
    return('phi_auto')
  }
  if(pgs_method == 'megaprs'){
    # Read in megaprs log file
    log <- readLines(paste0(outdir,'/reference/pgs_score_files/',pgs_method,'/',gwas,'/ref-',gwas,'.log'))
    log <- log[grepl('identified as the best with correlation', log)]
    pseudoval <- gsub(' .*','', gsub('Model ', '', log))
    return(paste0('ldak_Model', pseudoval))
  }
  if(pgs_method == 'lassosum'){
    # Read in megaprs log file
    log <- readLines(paste0(outdir,'/reference/pgs_score_files/',pgs_method,'/',gwas,'/ref-',gwas,'.log'))
    s_val <- gsub('.* ', '', log[grepl('^s = ', log)])
    lambda_val <- gsub('.* ', '', log[grepl('^lambda = ', log)])
    return(paste0('s', s_val, '_lambda', lambda_val))
  }

  # If pgs_method is external, return the only score
  if(pgs_method == 'external'){
    return('external')
  }
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

