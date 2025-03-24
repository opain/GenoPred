#!/usr/bin/Rscript

# Create the function for converting R2 into liability R2
h2l_R2 <- function(k, r2, p) {
  # K baseline disease risk
  # r2 from a linear regression model attributable to genomic profile risk score
  # P proportion of sample that are cases
  # calculates proportion of variance explained on the liability scale
  #from ABC at http://www.complextraitgenomics.com/software/
  #Lee SH, Goddard ME, Wray NR, Visscher PM. (2012) A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012 Apr;36(3):214-24.
  x= qnorm(1-k)
  z= dnorm(x)
  i=z/k
  C= k*(1-k)*k*(1-k)/(z^2*p*(1-p))
  theta= i*((p-k)/(1-k))*(i*((p-k)/(1-k))-x)
  h2l_R2 = C*r2 / (1 + C*theta*r2)
}

# Functions for reading in predictor file
read_predictor<-function(x, pred_miss, file_index = NULL, keep = NULL){
  # Read in predictor file
  tmp <- fread(x)
  
  # Create a column that combines FID and IID
  tmp <- combine_fid_iid(tmp)

  if(!is.null(keep)){
    # Restrict to keep
    tmp <- tmp[tmp$IID %in% keep,]
  }
  
  # Remove variables with > opt$pred_miss missing values
  tmp <- filter_columns_by_missing(tmp, threshold = opt$pred_miss)

  # Remove individuals with any missing data
  tmp <- tmp[complete.cases(tmp), ]

  # Update column names to avoid duplicate column names between predictor files
  names(tmp) <- gsub("[[:punct:]]", ".", names(tmp))
  names(tmp)[-1] <- paste0('PredFile', file_index, '.', names(tmp)[-1])

  log_add(log_file = log_file, message = paste0('Predictors file ', file_index, ' contains ', ncol(tmp) - 1, ' predictors with sufficient data.'))
  log_add(log_file = log_file, message = paste0('Predictors file ', file_index, ' contains ', nrow(tmp) ,' individuals with complete data for remaining predictors.'))

  return(tmp)
}

# Function for merging by IID
predictor_merger <- function(x, y) {
  return(merge(x, y, by = 'IID'))
}

# Function to combine FID and IID, into a variable called IID
combine_fid_iid <- function(data, fid_col = "FID", iid_col = "IID") {
  # Check if the specified columns exist in the data
  if (!all(c(fid_col, iid_col) %in% names(data))) {
    stop("Both FID and IID columns must be present in the data.")
  }

  # Combine FID and IID into a new IID column
  data[[iid_col]] <- paste0(data[[fid_col]], ':', data[[iid_col]])

  # Remove the FID column
  data[[fid_col]] <- NULL

  return(data)
}

# Remove variables with > opt$pred_miss missing values
filter_columns_by_missing <- function(data, threshold, first_col_keep = TRUE) {
  # Identify columns to keep based on missing value threshold
  col_keep <- c(first_col_keep, sapply(data[, -1, with = FALSE], function(col) {
    mean(!is.finite(col) | is.na(col)) < threshold
  }))

  # Filter the data by keeping only the columns that meet the criteria
  data_filtered <- data[, col_keep, with = FALSE]

  return(data_filtered)
}

# Set the grid search for the elastic net
# This is expanded from the default to allow more aggressive shrinkage
enet_grid <- expand.grid(
  alpha = seq(0, 1, length = 5),          # Explore alpha values: 0 (Ridge) to 1 (Lasso)
  lambda = 10^seq(-4, 1, length = 10)     # Explore lambda values: 0.0001 to 10
)

# Read in outcome file
read_outcome<-function(x, keep = NULL){
  outcome <- fread(x)
  outcome <- outcome[complete.cases(outcome), ]
  names(outcome)[3] <- 'outcome_var'

  # Create a column that combines FID and IID
  outcome <- combine_fid_iid(outcome)

  log_add(log_file = log_file, message = paste0('Outcome file contains ',nrow(outcome),' individuals with complete data.'))

  if(!is.null(keep)){
    ############
    # Extract individuals in the keep file
    ############

    # Read in keep file
    keep_file <- fread(keep)
    names(keep_file)[1:2] <- c('FID', 'IID')
    # Create a column that combines FID and IID
    keep_file <- combine_fid_iid(keep_file)
    # Extract keep individuals from the phenotypic data
    outcome <- outcome[(outcome$IID %in% keep_file$IID), ]

    log_add(log_file = log_file, message = paste0('Outcome file contains ',nrow(outcome),' individuals after extraction of individuals in ', keep,'.'))
  }

  return(outcome)
}

# Create a variable containing seeds for internal cross validation
fold_seeds<-function(n_fold){
  seeds <- vector(mode = "list", length = n_fold+1)
  for(i in 1:(n_fold)){
    seeds[[i]]<- sample.int(n=1000, 6)
  }
  seeds[[n_fold+1]]<-sample.int(n=1000, 1)

  return(seeds)
}

# Subset training and testing data
subset_train_test<-function(dat, train_ind, fold){

  dat_list<-list()
  dat_list$train<-list()
  dat_list$test<-list()

  # Subset data to training and testing sets
  dat_list$train$IID <- dat$IID[train_ind[[fold]]]
  dat_list$test$IID <- dat$IID[-train_ind[[fold]]]

  dat_list$train$y <- dat$outcome_var[train_ind[[fold]]]
  dat_list$test$y <- dat$outcome_var[-train_ind[[fold]]]

  dat_list$train$x <- dat[train_ind[[fold]], !(names(dat) %in% c('IID','outcome_var')), with=F]
  dat_list$test$x <- dat[-train_ind[[fold]], !(names(dat) %in% c('IID','outcome_var')), with=F]

  # If there is only one predictor, insert a column of 0s so elastic net function doesn't fail
  if(ncol(dat_list$train$x) == 1){
    dat_list$train$x<-data.table(cbind(0, dat_list$train$x))
    dat_list$test$x<-data.table(cbind(0, dat_list$test$x))
  }

  return(dat_list)
}

# Compare predicted and observed values
eval_pred <- function(obs, pred, family){
  mod <- summary(lm(scale(as.numeric(obs)) ~ scale(as.numeric(pred))))

  mod_sum <- data.table(
    R = coef(mod)[2, 1],
    SE = coef(mod)[2, 2],
    P = coef(mod)[2, 4]
  )

  if(family == 'binomial'){
    mod_sum <- data.table(
      mod_sum,
      R2l = h2l_R2(
        opt$outcome_pop_prev,
        coef(mod_sum)[2, 1] ^ 2,
        sum(obs == 'CASE') / length(obs)
      ),
      N = length(obs),
      Ncase = sum(obs == 'CASE'),
      Ncont = sum(obs == 'CONTROL')
    )
  } else {
    mod_sum <- data.table(
      mod_sum,
      R2o = coef(mod)[2, 1] ^ 2,
      N = length(obs)
    )
  }
  return(mod_sum)
}

#########
# Function for progress bar
#########

initialise_progress <- function(log_message, log_file){
  progress_file <- tempfile()
  saveRDS(0, progress_file)
  log_add(log_file = log_file, message = log_message)
  return(progress_file)
}

update_progress_file <- function(progress_file) {
  # Lock the file for writing (this ensures only one worker can modify it at a time)
  lockfile <- paste0(progress_file, ".lock")
  while (file.exists(lockfile)) {
    Sys.sleep(0.01)  # Wait if another worker is updating the file
  }

  # Create the lockfile
  file.create(lockfile)

  # Read current progress
  progress <- readRDS(progress_file)

  # Update progress_file
  progress <- progress + 1
  saveRDS(progress, progress_file)

  # Remove the lockfile
  file.remove(lockfile)

  return(progress)  # Return progress
}

update_log_file <- function(log_file, message) {
  # Read the current content of the log file
  log_content <- readLines(log_file)

  # Update the last line with the new message
  if (length(log_content) > 0) {
    log_content[length(log_content)] <- message
  } else {
    log_content <- message
  }

  # Write the updated log content back to the file
  writeLines(log_content, con = log_file)
}

export_final_model <- function(model, group, outdir){
  if ("glm" %in% class(model)) {
    model_coefficients <- as.matrix(coef(model))
    rownames(model_coefficients)[1]<-'intercept'
  }
  if ("glmnet" %in% class(model)) {
    model_coefficients <- coef(model, s = model$lambdaOpt)
    rownames(model_coefficients)[1]<-'intercept'
    model_coefficients <- as.matrix(model_coefficients[model_coefficients[,1] != 0, , drop=F])
  }

  write.table(
    model_coefficients,
    paste0(
      outdir,
      '/',
      group,
      '.final_model.txt'
    ),
    row.names = T,
    col.names = F,
    quote = F
  )
}
