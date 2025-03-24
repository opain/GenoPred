#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--outcome", action="store", default=NULL, type='character',
    help="File containing outcome data [required]"),
make_option("--predictors", action="store", default=NULL, type='character',
    help="File listing files containing predictors, with a groups column for model comparison [required]"),
make_option("--n_fold", action="store", default=10, type='numeric',
    help="Number of folds in for cross-validation [optional]"),
make_option("--n_core", action="store", default=1, type='numeric',
    help="Number of cores for parallel computing [optional]"),
make_option("--keep", action="store", default=NULL, type='character',
    help="File containing list of individuals to include in analysis [optional]"),
make_option("--outcome_pop_prev", action="store", default=NULL, type='numeric',
    help="Prevalence of outcome in the general population [optional]"),
make_option("--out", action="store", default=NULL, type='character',
    help="Prefix for output files [required]"),
make_option("--pred_miss", action="store", default=0.1, type='numeric',
    help="Proportion of missing values allowed in predictor [optional]"),
make_option("--export_models", action="store", default=T, type='logical',
    help="Export model coefficients [optional]"),
make_option("--seed", action="store", default=1, type='numeric',
    help="Set seed number [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load dependencies
suppressMessages(library(GenoUtils))
suppressMessages(library(data.table))
source('../functions/misc.R')
source_all('../functions')
suppressMessages(library(doMC))
suppressMessages(library(caret))
suppressMessages(library(pROC))
suppressMessages(library(verification))
suppressMessages(library(psych))
registerDoMC(opt$n_core)

# Check required inputs
if(is.null(opt$outcome)){
  stop('--outcome must be specified.\n')
}
if(is.null(opt$predictors)){
  stop('--predictors must be specified.\n')
}
if(is.null(opt$out)){
  stop('--out must be specified.\n')
}

# Create output directory
opt$output_dir <- paste0(dirname(opt$out), '/')
system(paste0('mkdir -p ', opt$output_dir))

# Create temp directory
tmp_dir <- tempdir()

# Create directory for final models to be saved
if(opt$export_models){
  system(paste0('mkdir -p ', opt$output_dir, '/final_models'))
}

# Initiate log file
log_file <- paste0(opt$out,'.log')
log_header(log_file = log_file, opt = opt, script = 'model_builder.R', start.time = start.time)

###########
# Read in the outcome data
###########

outcome<-read_outcome(x = opt$outcome, keep = opt$keep)

# Determine whether outcome is binary or continuous and format accordingly
if (length(unique(outcome$outcome_var)) > 2) {
  family <- 'gaussian'
}
if (length(unique(outcome$outcome_var)) == 2) {
  family <- 'binomial'
  outcome$outcome_var <- factor(outcome$outcome_var, labels = c('CONTROL', 'CASE'))
}

log_add(log_file = log_file, message = paste0('Phenotype is ', ifelse(family == 'binomial', 'binary', 'quantitative'),'.'))

###########
# Read in predictors
###########

predictors_file <- fread(opt$predictors)

if(nrow(predictors_file) > 1){
  predictors <- foreach(i = 1:nrow(predictors_file)) %dopar% {
    read_predictor(x = predictors_file$predictor[i], pred_miss = opt$pred_miss, file_index = i, keep = outcome$IID)
  }

  group_list <- do.call(rbind, lapply(1:nrow(predictors_file), function(predfile) {
    data.table(multi = predictors_file$multi[predfile], top1 = predictors_file$top1[predfile], predictor = names(predictors[[predfile]])[-1])
  }))

  predictors <- Reduce(function(x, y) merge(x, y, by = "IID"), predictors)

  log_add(log_file = log_file, message = paste0('After merging the ', nrow(predictors_file), ' predictors files, ', ncol(predictors)-1, ' predictors remain.'))
  log_add(log_file = log_file, message = paste0('After merging the ', nrow(predictors_file), ' predictors files, ', nrow(predictors), ' individuals remain.'))
} else {
  predictors <- read_predictor(x = predictors_file$predictor[1], pred_miss = opt$pred_miss)
  group_list <- data.table(multi = predictors_file$multi[1], top1 = predictors_file$top1[1], predictor = names(predictors)[-1])
}

# Remove predictors with zero variance
nz_var <- sapply(predictors[, -1, with = FALSE], function(col) var(col) != 0)
if(sum(!nz_var) > 1){
  log_add(log_file = log_file, message = paste0(sum(!nz_var), ' predictors have zero variance and will be excluded from downstream analyes.'))
}
if(all(!(nz_var))){
  stop('All predictors have zero variance.')
}
predictors <- predictors[, c(TRUE, nz_var), with = FALSE]
group_list <- group_list[group_list$predictor %in% names(predictors),]

###########
# Create list of groups for downstream comparison
###########

group_list$group <- paste0(group_list$multi,'-',group_list$top1)

# Remove identical predictors within each group
group_list_non_identical <- NULL
for(i in unique(group_list$group)){
  if(sum(group_list$group == i) > 1){
    ident <- group_list$predictor[group_list$group == i][
      duplicated(
        as.list(
          predictors[, group_list$predictor[group_list$group == i], with = F]))]

    group_list_non_identical <- rbind(
      group_list_non_identical,
      group_list[group_list$group == i & !(group_list$predictor %in% ident),]
    )

    if(length(ident) > 0){
      log_add(log_file = log_file, message = paste0(length(ident), ' duplicate predictors removed from group ', i))
    }
  } else {
    group_list_non_identical <- rbind(
      group_list_non_identical,
      group_list[group_list$group == i,]
    )
  }
}
group_list <- group_list_non_identical

# Calculate the number of predictors in each group
for(i in unique(group_list$multi)){
  group_list$n_multi[group_list$multi == i] <- sum(group_list$multi == i)
  for(j in unique(group_list$top1)){
    group_list$n_top1[group_list$multi == i & group_list$top1 == j] <- sum(group_list$multi == i & group_list$top1 == j)
  }
}

write.table(group_list[!duplicated(group_list$group), c('group','n_multi','n_top1'), with = F], paste0(opt$out,'.group_list.txt'), col.names=T, row.names=F, quote=F)
log_add(log_file = log_file, message = paste0('List of groups saved as ',opt$out,'.group_list.txt.'))

###########
# Merge the outcome and predictors
###########

outcome_predictors <- merge(outcome, predictors, by='IID')

rm(outcome, predictors)

log_add(log_file = log_file, message = paste0(nrow(outcome_predictors),' individuals have both phenotypic and predictor data.'))

# Report the size of the combined outcome and predictor data
log_add(log_file = log_file, message = paste0('Data to be carried foward is ',format(object.size(outcome_predictors), units='auto'),'.'))

############
# Prediction modelling
############

# We will derive top1 models, where the best predictor within each multi-top1 combo is evaluated using cross-validation
# We will then derive multi models, combining top1 predictors within each multi group, agian using cross-validation.

# Split the sample into opt$n_outer_fold folds
set.seed(opt$seed)
d <- sample(1:nrow(outcome_predictors))
train_ind <- createFolds(d, k = opt$n_fold, returnTrain=TRUE)

# Create objects to store outputs
indep_pred_list <- list()

####
# Generate predictions using single best predictor from each group, identifying the best predictor using training data, and then evaluating in the test data
####
# We will include the groups with a single predictor here as well for convenience, although no feature selection is required.

# Initialise progress log
log_message <- 'Generate predictions using top1 models for each group... '
progress_file <- initialise_progress(log_message = log_message, log_file = log_file)

top1_indep_pred <- foreach(i = 1:length(unique(group_list$group)), .combine = 'c') %dopar% {
  group_name<-paste0(unique(group_list$group)[i], '.top1')
  indep_pred <- NULL
  for(outer_val in 1:opt$n_fold){
    # Subset training and testing data
    cv_dat <- subset_train_test(dat = outcome_predictors, train_ind = train_ind, fold = outer_val)
    
    # Subset variables in group
    pred_name <- group_list$predictor[group_list$group == unique(group_list$group)[i]]
    cv_dat$train$x <- cv_dat$train$x[, pred_name, with = F]
    
    # Evaluate each predictor in training data
    # NOTE. Should we be using the RMSE to select the best predictor within a group.
    top1_res<-NULL
    for(pred_i in names(cv_dat$train$x)){
      res_pred_i <- cor(as.numeric(cv_dat$train$y), cv_dat$train$x[[pred_i]], use='p')
      top1_res <- rbind(
        top1_res,
        data.table(
          pred = pred_i,
          cor = res_pred_i)
      )
    }
    top_pred <- top1_res$pred[which.max(abs(top1_res$cor))]
    
    # Build model using best predictor
    train_tmp <- data.table(y = cv_dat$train$y, x = cv_dat$train$x[[top_pred]])
    train_mod <- glm(y ~ x, family=family, data=train_tmp)
    
    # Evaluate best performing predictor in test data
    test_tmp <- data.table(x = cv_dat$test$x[[top_pred]])
    indep_pred_i <- predict(object = train_mod, newdata = test_tmp, type = "response")
    indep_pred_i <- data.table(obs = cv_dat$test$y, pred = indep_pred_i)
    
    # Save test set predictions from each outer loop
    indep_pred <- rbind(indep_pred, indep_pred_i)
  }
  
  # Derive and export final model using all data
  if(opt$export_models){
    top1_res<-NULL
    for(pred_i in names(cv_dat$train$x)){
      res_pred_i <- cor(as.numeric(outcome_predictors$outcome_var), outcome_predictors[[pred_i]], use='p')
      top1_res <- rbind(
        top1_res,
        data.table(
          pred = pred_i,
          cor = res_pred_i)
      )
    }
    top_pred <- top1_res$pred[which.max(abs(top1_res$cor))]
    
    train_mod <- glm(as.formula(paste0('outcome_var ~ ', top_pred)), family=family, data=outcome_predictors)
    
    export_final_model(model = train_mod,
                       group = group_name,
                       outdir = paste0(opt$output_dir,
                                       '/final_models'))
  }
  
  # Update progress log
  progress <- update_progress_file(progress_file)
  update_log_file(log_file = log_file, message = paste0(log_message, floor(progress/length(unique(group_list$group))*100),'%'))
  
  # Output results
  setNames(list(indep_pred), group_name)
}

indep_pred_list <- c(indep_pred_list, top1_indep_pred)
update_log_file(log_file = log_file, message = paste0(log_message, 'Done!'))

####
# Generate predictions using model containing top1 predictores from each multi group
####

# Initialise progress log
log_message <- 'Generate predictions using multi models for each group... '
progress_file <- initialise_progress(log_message = log_message, log_file = log_file)

multi <- foreach(i = 1:length(unique(group_list$multi)), .combine = 'c') %dopar% {
  group_name<-paste0(unique(group_list$multi)[i], '.multi')
  indep_pred <- NULL
  for(outer_val in 1:opt$n_fold){
    # Subset training and testing data
    cv_dat <- subset_train_test(dat = outcome_predictors, train_ind = train_ind, fold = outer_val)
    
    top_pred_all <- NULL
    for(top1_group in unique(group_list$top1[group_list$multi == unique(group_list$multi)[i]])){
      # Subset variables in group
      pred_name <- group_list$predictor[
        group_list$top1 == top1_group & 
        group_list$multi == unique(group_list$multi)[i]]
      
      cv_dat_subset <- cv_dat
      cv_dat_subset$train$x <- cv_dat$train$x[, pred_name, with = F]
      
      # Evaluate each predictor in training data
      # NOTE. Should we be using the RMSE to select the best predictor within a group.
      top1_res<-NULL
      for(pred_i in names(cv_dat_subset$train$x)){
        res_pred_i <- cor(as.numeric(cv_dat_subset$train$y), cv_dat_subset$train$x[[pred_i]], use='p')
        top1_res <- rbind(
          top1_res,
          data.table(
            pred = pred_i,
            cor = res_pred_i)
        )
      }
      top_pred <- top1_res$pred[which.max(abs(top1_res$cor))]
      
      top_pred_all <- c(top_pred_all, top_pred)
    }
  
    # Build model using top1 predictors
    train_tmp <- data.table(y = cv_dat$train$y, cv_dat$train$x[,top_pred_all, with=F])
    train_mod <- glm(y ~ ., family=family, data=train_tmp)
    
    # Evaluate best performing predictor in test data
    test_tmp <- data.table(cv_dat$test$x[,top_pred_all, with=F])
    indep_pred_i <- predict(object = train_mod, newdata = test_tmp, type = "response")
    indep_pred_i <- data.table(obs = cv_dat$test$y, pred = indep_pred_i)
    
    # Save test set predictions from each outer loop
    indep_pred <- rbind(indep_pred, indep_pred_i)
  }
  
  # Derive and export final model using all data
  if(opt$export_models){
    top_pred_all <- NULL
    for(top1_group in unique(group_list$top1[group_list$multi == unique(group_list$multi)[i]])){
      # Subset variables in group
      pred_name <- group_list$predictor[
        group_list$top1 == top1_group & 
          group_list$multi == unique(group_list$multi)[i]]
      
      # Evaluate each predictor in training data
      # NOTE. Should we be using the RMSE to select the best predictor within a group.
      top1_res<-NULL
      for(pred_i in pred_name){
        res_pred_i <- cor(as.numeric(outcome_predictors$outcome_var), outcome_predictors[[pred_i]], use='p')
        top1_res <- rbind(
          top1_res,
          data.table(
            pred = pred_i,
            cor = res_pred_i)
        )
      }
      top_pred <- top1_res$pred[which.max(abs(top1_res$cor))]
      
      top_pred_all <- c(top_pred_all, top_pred)
    }
    
    # Build model using top1 predictors
    train_tmp <- data.table(y = outcome_predictors$outcome_var, outcome_predictors[,top_pred_all, with=F])
    train_mod <- glm(y ~ ., family=family, data=train_tmp)
    
    export_final_model(model = train_mod,
                       group = group_name,
                       outdir = paste0(opt$output_dir,
                                       '/final_models'))
  }
  
  # Update progress log
  progress <- update_progress_file(progress_file)
  update_log_file(log_file = log_file, message = paste0(log_message, floor(progress/length(unique(group_list$multi))*100),'%'))
  
  # Output results
  setNames(list(indep_pred), group_name)
}

indep_pred_list <- c(indep_pred_list, multi)
update_log_file(log_file = log_file, message = paste0(log_message, 'Done!'))

###################
# Evaluate all models
###################

# Initialise progress log
log_message <- 'Evaluating all models... '
progress_file <- initialise_progress(log_message = log_message, log_file = log_file)

pred_eval_all <- foreach(i = 1:length(indep_pred_list), .combine=rbind) %dopar% {
  # Update progress log
  progress <- update_progress_file(progress_file)
  update_log_file(log_file = log_file, message = paste0(log_message, floor(progress/length(indep_pred_list)*100),'%'))

  data.table(
      Group = names(indep_pred_list)[i],
      eval_pred(
        obs = indep_pred_list[[i]]$obs,
        pred = indep_pred_list[[i]]$pred,
        family = family
      )
    )
}

update_log_file(log_file = log_file, message = paste0(log_message, 'Done!'))

# Write out the results
write.table(pred_eval_all, paste0(opt$out,'.pred_eval.txt'), col.names=T, row.names=F, quote=F)
log_add(log_file = log_file, message = paste0('Model evaluation results saved as ',opt$out,'.pred_eval.txt.'))

###################
# Compare predictive utiliy of the different models
###################

if(length(pred_eval_all$Group) > 1){

  # Initialise progress log
  log_message <- 'Comparing model performance... '
  progress_file <- initialise_progress(log_message = log_message, log_file = log_file)

  # Identify number of pairwise comparisons
  models <- pred_eval_all$Group
  model_comps <- data.table(t(combn(models, 2)))

  comp_res_all <- foreach(i = 1:length(models), .combine=rbind) %dopar% {
    group1 <- models[i]
    comp_res_i <- NULL
    for(group2 in model_comps$V2[model_comps$V1 == group1]){

    	group1_r <- cor(scale(as.numeric(indep_pred_list[[group1]]$obs)), scale(indep_pred_list[[group1]]$pred))[1]
      group2_r <- cor(scale(as.numeric(indep_pred_list[[group2]]$obs)), scale(indep_pred_list[[group2]]$pred))[1]

      r_diff <- group1_r - group2_r

      group1_group2_r <- cor(scale(indep_pred_list[[group1]]$pred), scale(indep_pred_list[[group2]]$pred))

      r_diff_p <-
        paired.r(
          xy = group1_r,
          xz = group2_r,
          yz = group1_group2_r,
          n = length(scale(indep_pred_list[[group1]]$pred)),
          twotailed = T
        )$p[1]

      comp_res <- data.table(
        Model_1 = group1,
        Model_2 = group2,
        Model_1_R = group1_r,
        Model_2_R = group2_r,
        R_diff = r_diff,
        R_diff_pval = r_diff_p
      )

      comp_res_i <- rbind(comp_res_i, comp_res)
    }

    # Update progress log
    progress <- update_progress_file(progress_file)
    update_log_file(log_file = log_file, message = paste0(log_message, floor(progress/length(models)*100),'%'))

    comp_res_i
  }

  update_log_file(log_file = log_file, message = paste0(log_message, 'Done!'))

  # Write out the results
  write.table(comp_res_all, paste0(opt$out,'.pred_comp.txt'), col.names=T, row.names=F, quote=F)
  log_add(log_file = log_file, message = paste0('Model comparison results saved as ',opt$out,'.pred_comp.txt.'))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
log_add(log_file = log_file, message = paste0('Analysis finished at ',as.character(end.time)))
log_add(log_file = log_file, message = paste0('Analysis duration was ',as.character(round(time.taken,2)),attr(time.taken, 'units')))
