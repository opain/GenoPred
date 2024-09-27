#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--outcome", action="store", default=NULL, type='character',
    help="File containing outcome data [required]"),
make_option("--predictors", action="store", default=NULL, type='character',
    help="File listing files containing predictors, with a groups column for model comparison [required]"),
make_option("--n_outer_fold", action="store", default=10, type='numeric',
    help="Number of folds in for outer cross-validation [optional]"),
make_option("--n_inner_fold", action="store", default=10, type='numeric',
    help="Number of folds for inner cross-validation [optional]"),
make_option("--n_core", action="store", default=1, type='numeric',
    help="Number of cores for parallel computing [optional]"),
make_option("--keep", action="store", default=NULL, type='character',
    help="File containing list of individuals to include in analysis [optional]"),
make_option("--outcome_pop_prev", action="store", default=NULL, type='numeric',
    help="Prevalence of outcome in the general population [optional]"),
make_option("--out", action="store", default=NULL, type='character',
    help="Prefix for output files [required]"),
make_option("--assoc", action="store", default=T, type='logical',
    help="Perform association analysis between each predictor and outcome [optional]"),
make_option("--compare_predictors", action="store", default=F, type='logical',
    help="Option to assign each predictor to own group [optional]"),
make_option("--pred_miss", action="store", default=0.1, type='numeric',
    help="Proportion of missing values allowed in predictor [optional]"),
make_option("--top1", action="store", default=F, type='logical',
    help="Evaluate model using top predictor within each group [optional]"),
make_option("--all_model", action="store", default=T, type='logical',
    help="Evaluate model containing all predictors [optional]"),
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
suppressMessages(library(glmnet))
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
  system(paste0('mkdir ', opt$output_dir, '/final_models'))
}

# Initiate log file
log_file <- paste0(opt$out,'.log')
log_header(log_file = log_file, opt = opt, script = 'model_builder.R', start.time = start.time)

###########
# Read in predictors
###########
# We will create a table indicating which predictors belong to each group as well.

predictors_file <- fread(opt$predictors)

if(nrow(predictors_file) > 1){
  predictors <- foreach(i = 1:nrow(predictors_file)) %dopar% {
    read_predictor(x = predictors_file$predictor[i], pred_miss = opt$pred_miss, file_index = i)
  }

  group_list <- do.call(rbind, lapply(1:nrow(predictors_file), function(predfile) {
    data.table(group = predictors_file$group[predfile], predictor = names(predictors[[predfile]])[-1])
  }))

  predictors <- Reduce(function(x, y) merge(x, y, by = "IID"), predictors)

  log_add(log_file = log_file, message = paste0('After merging the ', nrow(predictors_file), ' predictors files, ', ncol(predictors)-1, ' predictors remain.'))
  log_add(log_file = log_file, message = paste0('After merging the ', nrow(predictors_file), ' predictors files, ', nrow(predictors), ' individuals remain.'))
} else {
  predictors <- read_predictor(x = predictors_file$predictor[1], pred_miss = opt$pred_miss)
  group_list <- data.table(group = predictors_file$group[1], predictor = names(predictors)[-1])
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

if(opt$compare_predictors){
  group_list <-
    rbind(group_list,
          data.table(
            group = paste0(group_list$group, '.', group_list$predictor),
            predictor = group_list$predictor
          ))
  log_add(log_file = log_file, message = 'Each predictor has been assigned to its own group.')
}


# Add a group containing all predictors
if(nrow(group_list) > 1 & opt$all_model){
  group_list <-
    rbind(group_list,
          data.table(
            group = 'all',
            predictor = unique(group_list$predictor)
          ))
  log_add(log_file = log_file, message = 'Group containing all predictors created.')
}

log_add(log_file = log_file, message = paste0(length(unique(group_list$group)), ' groups of predictors present.'))

# Calculate the number of predictors in each group
for(i in unique(group_list$group)){
  group_list$n[group_list$group == i] <- sum(group_list$group == i)
}

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
# Merge the outcome and predictors
###########

outcome_predictors <- merge(outcome, predictors, by='IID')

rm(outcome, predictors)

log_add(log_file = log_file, message = paste0(nrow(outcome_predictors),' individuals have both phenotypic and predictor data.'))

# Report the size of the combined outcome and predictor data
log_add(log_file = log_file, message = paste0('Data to be carried foward is ',format(object.size(outcome_predictors), units='auto'),'.'))

############
# Test association between outcome and each predictor
############

if(opt$assoc){
  log_add(log_file = log_file, message = 'Performing association analysis with each predictor...')

	outcome_predictors_y <- outcome_predictors$outcome_var
	outcome_predictors_x <- outcome_predictors[, -1:-2]

	assoc_res <- NULL
	assoc_res <- foreach(i = 1:ncol(outcome_predictors_x), .combine=rbind) %dopar% {

  	if(family == 'binomial'){
  	  mod <- glm(outcome_predictors_y ~ scale(outcome_predictors_x[[i]]), family=family)
  		obs_r2 <- cor(predict(mod), as.numeric(outcome_predictors_y))^2
  	  sum_mod <- summary(mod)
  	} else {
  	  mod <- glm(scale(outcome_predictors_y) ~ scale(outcome_predictors_x[[i]]), family = family)
  	  obs_r2 <- cor(predict(mod), as.numeric(outcome_predictors_y))^2
  	  sum_mod <- summary(mod)
  	}

	  data.table(
	    Group = group_list$group[group_list$predictor == names(outcome_predictors_x)[i]][1],
	    Predictor = names(outcome_predictors_x)[i],
	    BETA = coef(sum_mod)[2, 1],
	    SE = coef(sum_mod)[2, 2],
	    P = coef(sum_mod)[2, 4],
	    Obs_R2 = obs_r2,
  		N = length(outcome_predictors_y)
	  )

	}

	if(family == 'binomial'){
  	assoc_res$N_case <- sum(outcome_predictors_y == 'CASE')
  	assoc_res$N_control <- sum(outcome_predictors_y == 'CONTROL')
  	assoc_res$Liab_R2 <- h2l_R2(
  	  opt$outcome_pop_prev,
  	  assoc_res$obs_r2,
  	  sum(outcome_predictors_y == 'CASE') / length(outcome_predictors_y)
  	)
	}

	# Write out the results
	write.table(assoc_res, paste0(opt$out,'.assoc.txt'), col.names=T, row.names=F, quote=F)
	log_add(log_file = log_file, message = paste0('Predictor association results saved as ',opt$out,'.assoc.txt.'))
}

############
# Prediction modelling
############
# We will use elastic net model when groups contain more than one predictor, and a glm when the group contains only 1 predictor.
# Elastic net models will be derived and evaluated using nested cross validation, but glm will be evaluated using standard cross validation.
# We will store the predicted and observed values from the outer loops for each model, and then evaluate the models at the end

# Split the sample into opt$n_outer_fold folds
set.seed(opt$seed)
d <- sample(1:nrow(outcome_predictors))
train_ind <- createFolds(d, k = opt$n_outer_fold, returnTrain=TRUE)

# Set seeds for internal loop of nested CV for elastic net
seeds <- fold_seeds(opt$n_outer_fold)

# Create objects to store outputs
indep_pred_list <- list()

####
# Generate predictions using single best predictor from each group, identifying the best predictor using training data, and then evaluating in the test data
####

if(opt$top1){
  # Only run for groups containing more than 1 predictor
  top1_groups <- unique(group_list$group[group_list$n > 1])

  # Initialise progress log
  log_message <- 'Generate predictions using top1 models for each group... '
  progress_file <- initialise_progress(log_message = log_message, log_file = log_file)

  top1_indep_pred <- foreach(i = 1:length(top1_groups), .combine = 'c') %dopar% {
    group_name<-paste0(top1_groups[i], '.top1')
    indep_pred <- NULL
    for(outer_val in 1:opt$n_outer_fold){
      # Subset training and testing data
      cv_dat <- subset_train_test(dat = outcome_predictors, train_ind = train_ind, fold = outer_val)

      # Subset variables in group
      pred_name <- group_list$predictor[group_list$group == top1_groups[i]]
      cv_dat$train$x <- cv_dat$train$x[, pred_name, with = F]

      # Evaluate each predictor in training data
      # NOTE. Should we be using the RMSE to select the best predictor within a group.
      top1_res<-NULL
      for(pred_i in names(cv_dat$train$x)){
        res_pred_i <- cor(cv_dat$train$y, cv_dat$train$x[[pred_i]], use='p')
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
        res_pred_i <- cor(outcome_predictors$outcome_var, outcome_predictors[[pred_i]], use='p')
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
    update_log_file(log_file = log_file, message = paste0(log_message, floor(progress/length(top1_groups)*100),'%'))

    # Output results
    setNames(list(indep_pred), group_name)
  }

  indep_pred_list <- c(indep_pred_list, top1_indep_pred)
  update_log_file(log_file = log_file, message = paste0(log_message, 'Done!'))
}

###
# Generate predictions using single predictor groups
###

if(any(group_list$n == 1)){
  single_groups <- unique(group_list$group[group_list$n == 1])

  # Initialise progress log
  log_message <- 'Generate predictions using models containing single predictor... '
  progress_file <- initialise_progress(log_message = log_message, log_file = log_file)

  single_indep_pred <- foreach(i = 1:length(single_groups), .combine=c) %dopar% {
    indep_pred <- NULL
    for(outer_val in 1:opt$n_outer_fold){
      # Subset training and testing data
      cv_dat <- subset_train_test(dat = outcome_predictors, train_ind = train_ind, fold = outer_val)

      # Subset variables in group
      pred_name <- group_list$predictor[group_list$group == single_groups[i]]
      cv_dat$train$x <- cv_dat$train$x[[pred_name]]

      # Build model using predictor
      train_tmp<-data.table(y = cv_dat$train$y, x = cv_dat$train$x)
      train_mod<-glm(y ~ x, family=family, data=train_tmp)

      # Evaluate best performing predictor in test data
      test_tmp <- data.table(x = cv_dat$test$x[[pred_name]])
      indep_pred_i <- predict(object = train_mod, newdata = test_tmp, type = "response")
      indep_pred_i <- data.table(obs = cv_dat$test$y, pred = indep_pred_i)

      # Save test set predictions from each outer loop
      indep_pred <- rbind(indep_pred, indep_pred_i)
    }

    # Derive and export final model using all data
    if(opt$export_models){
      train_mod <- glm(as.formula(paste0('outcome_var ~ ', pred_name)), family=family, data=outcome_predictors)

      export_final_model(model = train_mod,
                         group = single_groups[i],
                         outdir = paste0(opt$output_dir,
                                         '/final_models'))
    }

    # Update progress log
    progress <- update_progress_file(progress_file)
    update_log_file(log_file = log_file, message = paste0(log_message, floor(progress/length(single_groups)*100),'%'))

    # Output results
    setNames(list(indep_pred), single_groups[i])
  }

  indep_pred_list <- c(indep_pred_list, single_indep_pred)
  update_log_file(log_file = log_file, message = paste0(log_message, 'Done!'))
}

###
# Generate predictions using elastic net model containing all predictors within each group
###
# Only apply to groups with more than one predictor
if(any(group_list$n > 1)){

  # Initialise progress log
  log_message <- 'Generate predictions using elastic net models for each group... '
  progress_file <- initialise_progress(log_message = log_message, log_file = log_file)

  group_enet <- unique(group_list$group[group_list$n > 1])
  enet_indep_pred <- list()
  for(i in 1:length(group_enet)){
    indep_pred <- NULL
    for(outer_val in 1:opt$n_outer_fold){
      # Subset training and testing data
      cv_dat <- subset_train_test(dat = outcome_predictors, train_ind = train_ind, fold = outer_val)

      # Subset variables in group
      pred_name <- group_list$predictor[group_list$group == group_enet[i]]
      cv_dat$train$x <- cv_dat$train$x[, pred_name, with = F]
      cv_dat$test$x <- cv_dat$test$x[, pred_name, with = F]

      # Train elastic net
      model <-
        train(
          y = cv_dat$train$y,
          x = cv_dat$train$x,
          trControl = trainControl(
            method = "cv",
            seeds = seeds,
            number = opt$n_inner_fold
          ),
          method = "glmnet",
          family = family,
          tuneGrid = enet_grid
        )

      indep_pred_i <- as.numeric(predict(object = model$finalModel, newx = data.matrix(cv_dat$test$x), type = "response", s = model$finalModel$lambdaOpt))
      indep_pred_i <- data.table(obs = cv_dat$test$y, pred = indep_pred_i)

      # Save test set predictions from each outer loop
      indep_pred <- rbind(indep_pred, indep_pred_i)
    }

    # Derive and export final model using all data
    if(opt$export_models){
      model <-
        train(
          y = outcome_predictors$outcome_var,
          x = outcome_predictors[, pred_name, with=F],
          trControl = trainControl(
            method = "cv",
            seeds = seeds,
            number = opt$n_inner_fold
          ),
          method = "glmnet",
          family = family,
          tuneGrid = enet_grid
        )

      export_final_model(model = model$finalModel,
                         group = group_enet[i],
                         outdir = paste0(opt$output_dir,
                                         '/final_models'))
    }

    # Update progress log
    progress <- update_progress_file(progress_file)
    update_log_file(log_file = log_file, message = paste0(log_message, floor(progress/length(group_enet)*100),'%'))

    # Output results
    enet_indep_pred[[group_enet[i]]] <- indep_pred
  }

  indep_pred_list <- c(indep_pred_list, enet_indep_pred)
  update_log_file(log_file = log_file, message = paste0(log_message, 'Done!'))
}

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
log_add(log_file = log_file, message = paste0('Model evaluation results saved as ',opt$out,'.pred_comp.txt.'))

end.time <- Sys.time()
time.taken <- end.time - start.time
log_add(log_file = log_file, message = paste0('Analysis finished at ',as.character(end.time)))
log_add(log_file = log_file, message = paste0('Analysis duration was ',as.character(round(time.taken,2)),attr(time.taken, 'units')))
