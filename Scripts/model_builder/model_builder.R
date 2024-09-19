#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--pheno", action="store", default=NULL, type='character',
		help="File containing phenotypic data [required]"),
make_option("--predictors", action="store", default=NULL, type='character',
		help="File listing files containing predictors, with a groups column for model comparison [required]"),
make_option("--n_outer_fold", action="store", default=5, type='numeric',
    help="Number of folds in for outer cross-validation [optional]"),
make_option("--n_inner_fold", action="store", default=10, type='numeric',
    help="Number of folds for inner cross-validation [optional]"),
make_option("--n_core", action="store", default=1, type='numeric',
		help="Number of cores for parallel computing [optional]"),
make_option("--keep", action="store", default=NULL, type='character',
		help="File containing list of individuals to include in analysis [optional]"),
make_option("--outcome_pop_prev", action="store", default=NULL, type='numeric',
		help="Prevalence of outcome in the general population [optional]"),
make_option("--output", action="store", default=NULL, type='character',
		help="Prefix for output files [required]"),
make_option("--save_group_model", action="store", default=F, type='logical',
		help="Save group models for external validation [optional]"),
make_option("--assoc", action="store", default=T, type='logical',
		help="Perform association analysis between each predictor and outcome [optional]"),
make_option("--compare_predictors", action="store", default=F, type='logical',
    help="Option to assign each predictor to own group [optional]"),
make_option("--interaction", action="store", default=F, type='logical',
    help="Option to include interaction terms between predictors in each group [optional]"),
make_option("--pred_miss", action="store", default=0.1, type='numeric',
    help="Proportion of missing values allowed in predictor [optional]"),
make_option("--top1", action="store", default=F, type='numeric',
    help="Evaluate model using top predictor within each group [optional]"),
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
if(is.null(opt$pheno)){
  stop('--pheno must be specified.\n')
}
if(is.null(opt$predictors)){
  stop('--predictors must be specified.\n')
}
if(is.null(opt$output)){
  stop('--out must be specified.\n')
}

# Create output directory
opt$output_dir <- paste0(dirname(opt$output), '/')
system(paste0('mkdir -p ', opt$output_dir))

# Create temp directory
tmp_dir <- tempdir()

# Initiate log file
log_file <- paste0(opt$output,'.log')
log_header(log_file = log_file, opt = opt, script = 'model_builder.R', start.time = start.time)

##################

# Read in the predictors file
predictors_list <- fread(opt$predictors)
names(predictors_list)[1] <- 'predictor'

# Determine if there is a group column, and if there is more than 1 group.
if (opt$compare_predictors == F) {
  if (!is.null(predictors_list$group)) {
    if (length(unique(predictors_list$group)) != 1) {
      opt$group_info <- T
      opt$model_comp <- T
      predictors_list$group <-
        gsub("[[:punct:]]", ".", predictors_list$group)
      log_add(log_file = log_file, message = 'Predictors file contains group information so model comparisons will be performed.')
    } else {
      opt$group_info <- T
      opt$model_comp <- F
      log_add(log_file = log_file, message = 'Predictors file does not contain group information so model comparisons will not be performed.')
    }
  } else {
    opt$group_info <- F
    opt$model_comp <- F
    log_add(log_file = log_file, message = 'Predictors file does not contain group information so model comparisons will not be performed.')
  }
} else {
  if (!is.null(predictors_list$group)) {
    opt$group_info <- T
    if (length(unique(predictors_list$group)) != 1) {
      predictors_list$group <-
        gsub("[[:punct:]]", ".", predictors_list$group)
      log_add(log_file = log_file, message = 'Predictors file contains group information so model comparisons will be performed.')
    } else {
      log_add(log_file = log_file, message = 'Predictors file does not contain group information.')
    }
  } else {
    opt$group_info <- F
    log_add(log_file = log_file, message = 'Predictors file does not contain group information.')
  }
  log_add(log_file = log_file, message = 'Each predictor will be assigned to a seperate group and model comparisons will be performed.')
  opt$model_comp <- T
}

###########
# Read in the phenotypic data
###########

Outcome <- fread(opt$pheno)
Outcome <- Outcome[complete.cases(Outcome), ]
names(Outcome)[3] <- 'Outcome_var'

# Create a column that combines FID and IID
Outcome <- combine_fid_iid(Outcome)

log_add(log_file = log_file, message = paste0('Phenotype file contains ',nrow(Outcome),' individuals with complete data.'))

# Determine whether outcome is binary or continuous and format accordingly
if (length(unique(Outcome$Outcome_var)) > 2) {
  opt$family <- 'gaussian'
}
if (length(unique(Outcome$Outcome_var)) == 2) {
  opt$family <- 'binomial'
  Outcome$Outcome_var <- factor(Outcome$Outcome_var, labels = c('CONTROL', 'CASE'))
}

log_add(log_file = log_file, message = paste0('Phenotype is ', ifelse(opt$family == 'binomial', 'binary', 'quantitative'),'.'))

if(!is.null(opt$keep)){
	############
	# Extract individuals in the keep file
	############

  # Read in keep file
  keep_file <- fread(opt$keep)
  names(keep_file)[1:2] <- c('FID', 'IID')
  # Create a column that combines FID and IID
  keep_file <- combine_fid_iid(keep_file)
  # Extract keep individuals from the phenotypic data
  Outcome <- Outcome[(Outcome$IID %in% keep_file$IID), ]

	log_add(log_file = log_file, message = paste0('Phenotype file contains ',nrow(Outcome),' individuals after extraction of individuals in ',opt$keep,'.'))
}

############
# Read in the predictor variables
############

if(nrow(predictors_list) > 1) {
  Predictors<-foreach(k = 1:dim(predictors_list)[1], .combine = predictor_merger) %dopar% {
    # Read in predictor file k
    tmp <- fread(predictors_list$predictor[k])

    # Create a column that combines FID and IID
    tmp <- combine_fid_iid(tmp)

    # Remove variables with > opt$pred_miss missing values
    tmp <- filter_columns_by_missing(tmp, threshold = opt$pred_miss)

    # Remove individuals with any missing data
    tmp <- tmp[complete.cases(tmp), ]

    # Update column names to avoid duplicate column names between predictor files
    names(tmp)[-1] <- paste0('PredFile', k, '.', names(tmp)[-1])

    if(opt$model_comp == T){
      # Add group name to each predictor
      names(tmp)[-1]<-paste0('Group_',predictors_list$group[k],'.',names(tmp)[-1])
    }

    log_add(log_file = log_file, message = paste0('Predictors file ',k,' contains ',ncol(tmp)-1,' predictors with sufficient data.'))
    log_add(log_file = log_file, message = paste0('Predictors file ',k,' contains ',nrow(tmp),' individuals with complete data for remaining predictors.'))

    tmp
  }

  Predictors<-data.table(Predictors)

  log_add(log_file = log_file, message = paste0('After merging the ', length(predictors_list), ' predictors files,', ncol(Predictors)-1, ' predictors remain.'))
  log_add(log_file = log_file, message = paste0('After merging the ', length(predictors_list), ' predictors files,', nrow(Predictors), ' individuals remain.'))

} else {

	# Read in the predictor variables
	Predictors<-fread(predictors_list$predictor[1])

	# Create a column that combines FID and IID
	Predictors <- combine_fid_iid(Predictors)

	# Remove variables with > opt$pred_miss missing values
	Predictors <- filter_columns_by_missing(Predictors, threshold = opt$pred_miss)

	# Remove individuals with any missing data
	Predictors<-Predictors[complete.cases(Predictors),]

	if(opt$compare_predictors == T){
		# Create the object predictors_list
		predictors_list <- data.table(predictor=names(Predictors)[-1], group=names(Predictors)[-1])
		predictors_list$group <- gsub("[[:punct:]]", ".",predictors_list$group)

		# Add group name to each predictor
		names(Predictors)[-1]<-paste0('Group_',names(Predictors)[-1],'.',names(Predictors)[-1])
	}

	log_add(log_file = log_file, message = paste0('Predictors file contains ', ncol(Predictors)-1, ' predictors with sufficient data.'))
	log_add(log_file = log_file, message = paste0('Predictors file contains ', dim(Predictors)[1], ' individuals with complete data for remaining predictors.'))

}

predictors_list_new<-predictors_list

if(opt$compare_predictors == T){
  if(opt$group_info == T){
    # Create the object predictors_list
    predictors_list_new <-
      rbind(predictors_list_new,
            data.table(
              predictor = names(Predictors)[-1],
              group = names(Predictors)[-1]
            ))
    predictors_list_new$group <- gsub("[[:punct:]]", ".",predictors_list_new$group)
    predictors_list_new$group <- gsub('Group.','', predictors_list_new$group)
  }
}

###########
# Merge the phenotype and predictor variables
###########

Outcome_Predictors <- merge(Outcome, Predictors, by='IID')

rm(Outcome,Predictors)

log_add(log_file = log_file, message = paste0(nrow(Outcome_Predictors),' individuals have both phenotypic and predictor data.'))

# Report the size of the combined outcome and predictor data
log_add(log_file = log_file, message = paste0('Data to be carried foward is ',format(object.size(Outcome_Predictors), units='auto'),'.'))

if(opt$assoc == T){

	############
	# Test association between Outcome and each variable in Predictors
	############

  log_add(log_file = log_file, message = 'Performing association analysis with each predictor...')

	Outcome_Predictors_y <- Outcome_Predictors$Outcome_var
	Outcome_Predictors_x <- Outcome_Predictors[, -1:-2]

	Assoc_res<-NULL

	if(opt$family == 'binomial'){
		Assoc_res<-foreach(i=1:ncol(Outcome_Predictors_x), .combine=rbind) %dopar% {
			if(var(Outcome_Predictors_x[[i]]) == 0){
				Assoc_res_temp<-data.frame(	Predictor = names(Outcome_Predictors_x)[i],
																		BETA = NA,
																		SE = NA,
																		P = NA,
																		Obs_R2 = NA)

				Assoc_res_temp
			} else {
			  mod <- glm(Outcome_Predictors_y ~ scale(Outcome_Predictors_x[[i]]), family=opt$family)
				obs_r2 <- cor(predict(mod), as.numeric(Outcome_Predictors_y))^2
			  sum_mod <- summary(mod)
				prob <- predict(mod,type=c("response"))
				Assoc_res_temp <-
				  data.table(
				    Predictor = names(Outcome_Predictors_x)[i],
				    BETA = coef(sum_mod)[2, 1],
				    SE = coef(sum_mod)[2, 2],
				    P = coef(sum_mod)[2, 4],
				    Obs_R2 = obs_r2
				  )
			  Assoc_res_temp
			}
		}
		# Convert Nagelkerke R2 to liability scale
		Assoc_res$Liab_R2 <-
		  h2l_R2(
		    opt$outcome_pop_prev,
		    Assoc_res$Obs_R2,
		    sum(Outcome_Predictors_y == 'CASE') / length(Outcome_Predictors_y)
		  )
	} else {
		Assoc_res <- foreach(i = 1:ncol(Outcome_Predictors_x), .combine=rbind) %dopar% {
			if(var(Outcome_Predictors_x[[i]]) == 0){
				Assoc_res_temp<-data.frame(	Predictor=names(Outcome_Predictors_x)[i],
																		BETA=NA,
																		SE=NA,
																		P=NA,
																		Obs_R2=NA)

				Assoc_res_temp
			} else {
				mod <- glm(scale(Outcome_Predictors_y) ~ scale(Outcome_Predictors_x[[i]]), family = opt$family)
				sum_mod <- summary(mod)
				Assoc_res_temp <-
				  data.table(
				    Predictor = names(Outcome_Predictors_x)[i],
				    BETA = coef(sum_mod)[2, 1],
				    SE = coef(sum_mod)[2, 2],
				    P = coef(sum_mod)[2, 4],
				    Obs_R2 = coef(sum_mod)[2, 1] ^ 2
				  )
				Assoc_res_temp
			}
		}
	}

	# Write out the results
	write.table(Assoc_res, paste0(opt$output,'.assoc.txt'), col.names=T, row.names=F, quote=F)
	log_add(log_file = log_file, message = paste0('Predictor association results saved as ',opt$output,'.assoc.txt.'))
}

if(opt$model_comp == T){
  ############
  # Build and evaluate models using predictors together
  ############
  # We will use nested cross validation to evaluate models

  # Split the sample into n=opt$internal_validation_prop equal parts
  set.seed(opt$seed)
  nr <- nrow(Outcome_Predictors)
  d <- sample(1:nr)

  train.ext = createFolds(d, k = opt$n_outer_fold, returnTrain = TRUE)
  test.ext = lapply(train.ext, function(x) (1:nr)[-x])

  # Create a variable containing seeds for internal cross validation
  seeds <- vector(mode = "list", length = opt$n_inner_fold+1)
  for(i in 1:(opt$n_inner_fold)){
    seeds[[i]] <- sample.int(n=1000, 10)
  }
  seeds[[opt$n_inner_fold + 1]] <- sample.int(n = 1000, 1)

  # Create a model for each group, and using all predictors
  predictors_list_new <-
    rbind(predictors_list_new,
          data.frame(predictor = 'All',
                     group = '.'))

  log_add(log_file = log_file, message = 'Initiating nested cross-validation...')

  indep_pred<-list()
  Prediction_summary_all<-NULL

  if(opt$model_comp == T){
    if(opt$top1){
      ##############
      # Evaluate top1 model
      ##############
      for(group in unique(predictors_list_new$group[predictors_list_new$group != '.'])){
        group_name<-paste0(group, '_top1')
        for(outer_val in 1:opt$n_outer_fold){
          Outcome_Predictors_train_ind <- d[train.ext[[outer_val]]]

          Outcome_Predictors_train <- Outcome_Predictors[Outcome_Predictors_train_ind,]
          Outcome_Predictors_test <- Outcome_Predictors[-Outcome_Predictors_train_ind,]

          Outcome_Predictors_train_y <- Outcome_Predictors_train$Outcome_var
          Outcome_Predictors_train_x <- Outcome_Predictors_train[, -1:-2]

          Outcome_Predictors_test_y <- Outcome_Predictors_test$Outcome_var
          Outcome_Predictors_test_x <- Outcome_Predictors_test[, -1:-2]

          Outcome_Predictors_train_x_group <-
            Outcome_Predictors_train_x[,grepl(
              paste0('Group_', group, '\\.', '|', 'Group_', group, '$'),
              names(Outcome_Predictors_train_x)
            ), with=F]

          # Evaluate each predictor in training data
          top1_res<-NULL
          for(pred_i in names(Outcome_Predictors_train_x_group)){
            res_pred_i<-cor(Outcome_Predictors_train_y, Outcome_Predictors_train_x_group[[pred_i]], use='p')
            top1_res<-rbind(
              top1_res,
              data.frame(
                pred=pred_i,
                cor=res_pred_i)
            )
          }

          # Build model using best predictor
          train_tmp<-data.frame(
            y=Outcome_Predictors_train_y,
            x=Outcome_Predictors_train_x_group[[top1_res$pred[which.max(abs(top1_res$cor))]]])
          train_mod<-glm(y ~ x, family=opt$family, data=train_tmp)

          # Evaluate best performing predictor in test data
          test_tmp<-data.frame(
            x=Outcome_Predictors_test_x[[top1_res$pred[which.max(abs(top1_res$cor))]]])

          Indep_Pred<-predict(object = train_mod, newdata = test_tmp, type = "response")

          rm(train_tmp)
          rm(test_tmp)
          Indep_Pred_tab<-data.frame(obs=Outcome_Predictors_test_y,
                                     pred=Indep_Pred)

          print(top1_res$pred[which.max(abs(top1_res$cor))])
          print(cor(Indep_Pred_tab))

          # Save test set predictions from each outer loop
          indep_pred[[group_name]]<-rbind(indep_pred[[group_name]], Indep_Pred_tab)

        }
        Indep_mod<-summary(lm(scale(as.numeric(indep_pred[[group_name]]$obs)) ~ scale(indep_pred[[group_name]]$pred)))

        if(opt$family=='binomial'){
          Prediction_summary_all<-rbind(Prediction_summary_all, data.frame(	Model=paste0(group_name,'_group'),
                                                                            R=coef(Indep_mod)[2,1],
                                                                            SE=coef(Indep_mod)[2,2],
                                                                            P=coef(Indep_mod)[2,4],
                                                                            R2l=h2l_R2(opt$outcome_pop_prev, coef(Indep_mod)[2,1]^2, sum(Outcome_Predictors_train_y == 'CASE')/length(Outcome_Predictors_train_y)),
                                                                            N=dim(indep_pred[[group_name]])[1],
                                                                            Ncase=sum(Outcome_Predictors_train_y == 'CASE'),
                                                                            Ncont=sum(Outcome_Predictors_train_y == 'CONTROL')))
        } else {
          Prediction_summary_all<-rbind(Prediction_summary_all, data.frame(	Model=paste0(group_name,'_group'),
                                                                            R=coef(Indep_mod)[2,1],
                                                                            SE=coef(Indep_mod)[2,2],
                                                                            P=coef(Indep_mod)[2,4],
                                                                            R2o=coef(Indep_mod)[2,1]^2,
                                                                            N=dim(indep_pred[[group_name]])[1]))
        }
      }
    }

    # Build glmnet using each group of predictors at a time
    for(group in unique(predictors_list_new$group)){
      for(outer_val in 1:opt$n_outer_fold){
        print(outer_val)

        Outcome_Predictors_train_ind <- d[train.ext[[outer_val]]]

        Outcome_Predictors_train <- Outcome_Predictors[Outcome_Predictors_train_ind,]
        Outcome_Predictors_test <- Outcome_Predictors[-Outcome_Predictors_train_ind,]

        Outcome_Predictors_train_y <- Outcome_Predictors_train$Outcome_var
        Outcome_Predictors_train_x <- Outcome_Predictors_train[, -1:-2]

        Outcome_Predictors_test_y <- Outcome_Predictors_test$Outcome_var
        Outcome_Predictors_test_x <- Outcome_Predictors_test[, -1:-2]

        # Subset predictor in the group
        print(group)
        Outcome_Predictors_train_x_group <-
          Outcome_Predictors_train_x[,grepl(
            paste0('Group_', group, '\\.', '|', 'Group_', group, '$'),
            names(Outcome_Predictors_train_x)
          ), with=F]
        Outcome_Predictors_test_x_group <-
          Outcome_Predictors_test_x[,grepl(
            paste0('Group_', group, '\\.', '|', 'Group_', group, '$'),
            names(Outcome_Predictors_test_x)
          ), with=F]

    		# Rename group to 'all' if == '.'
    		if (group == '.') {
    		  group_name <- 'All'
    		  Outcome_Predictors_train_x_group <- Outcome_Predictors_train_x
    		  Outcome_Predictors_test_x_group <- Outcome_Predictors_test_x
    		} else {
    		  group_name <- group
    		}

    		# If there is only one predictor, use glm
    		if(dim(Outcome_Predictors_train_x_group)[2] > 1){
    				model<- train(y=Outcome_Predictors_train_y, x=Outcome_Predictors_train_x_group, trControl=trainControl(method="cv", seeds=seeds, number=opt$n_inner_fold, classProbs=T, savePredictions = 'final'), method="glmnet", family=opt$family)
    		} else {
    				model<- train(y=Outcome_Predictors_train_y, x=cbind(0,Outcome_Predictors_train_x_group), trControl=trainControl(method="cv", seeds=seeds, number=opt$n_inner_fold, classProbs=T, savePredictions = 'final'), method="glm", family=opt$family)
    		}

  		  if(dim(Outcome_Predictors_train_x_group)[2] > 1){
  				Indep_Pred<-as.numeric(predict(object = model$finalModel, newx = data.matrix(Outcome_Predictors_test_x_group), type = "response", s = model$finalModel$lambdaOpt))
  			} else {
  					tmp<-data.frame(cbind(0,Outcome_Predictors_test_x_group))
  					names(tmp)[1]<-'0'
  					Indep_Pred<-predict(object = model$finalModel, newdata = tmp, type = "response")
  					rm(tmp)
  			}

  		  Indep_Pred_tab<-data.frame(obs=Outcome_Predictors_test_y,
  		                             pred=Indep_Pred)

  			# Save test set predictions from each outer loop
  			indep_pred[[group_name]]<-rbind(indep_pred[[group_name]], Indep_Pred_tab)
      }

      Indep_mod<-summary(lm(scale(as.numeric(indep_pred[[group_name]]$obs)) ~ scale(indep_pred[[group_name]]$pred)))

      if(opt$family=='binomial'){
        Prediction_summary_all<-rbind(Prediction_summary_all, data.frame(	Model=paste0(group_name,'_group'),
                                                                          R=coef(Indep_mod)[2,1],
                                                                          SE=coef(Indep_mod)[2,2],
                                                                          P=coef(Indep_mod)[2,4],
                                                                          R2l=h2l_R2(opt$outcome_pop_prev, coef(Indep_mod)[2,1]^2, sum(Outcome_Predictors_train_y == 'CASE')/length(Outcome_Predictors_train_y)),
                                                                          N=dim(indep_pred[[group_name]])[1],
                                                                          Ncase=sum(Outcome_Predictors_train_y == 'CASE'),
                                                                          Ncont=sum(Outcome_Predictors_train_y == 'CONTROL')))
      } else {
        Prediction_summary_all<-rbind(Prediction_summary_all, data.frame(	Model=paste0(group_name,'_group'),
                                                                          R=coef(Indep_mod)[2,1],
                                                                          SE=coef(Indep_mod)[2,2],
                                                                          P=coef(Indep_mod)[2,4],
                                                                          R2o=coef(Indep_mod)[2,1]^2,
                                                                          N=dim(indep_pred[[group_name]])[1]))
      }
    }

    # Write out the results
    write.table(Prediction_summary_all, paste0(opt$output,'.pred_eval.txt'), col.names=T, row.names=F, quote=F)

    log_add(log_file = log_file, message = paste0('Model evaluation results saved as ',opt$output,'.pred_eval.txt.'))

  	###################
  	# Compare predictive utiliy of the different models
  	###################
    if(opt$top1){
      predictors_list_new_top1<-predictors_list_new[-nrow(predictors_list_new),]
      predictors_list_new_top1$group<-paste0(predictors_list_new_top1$group, '_top1')
      predictors_list_new<-rbind(predictors_list_new, predictors_list_new_top1)
    }
  	comp_res_all<-NULL
  	for(group1 in unique(predictors_list_new$group)){
    	for(group2 in unique(predictors_list_new$group)){
    	    if(group1 == '.'){
    	      group1<-'All'
    	    }
      	  if(group2 == '.'){
    	      group2<-'All'
      	  }
    	    if(group1 == group2){
      			group1_r<-cor(scale(as.numeric(indep_pred[[group1]]$obs)), scale(indep_pred[[group1]]$pred))[1]
    			  group2_r<-cor(scale(as.numeric(indep_pred[[group2]]$obs)), scale(indep_pred[[group2]]$pred))[1]

    	      comp_res<-data.frame(Model_1=group1,
    	                           Model_2=group2,
    	                           Model_1_R=group1_r,
    	                           Model_2_R=group2_r,
    	                           R_diff=0,
    	                           R_diff_pval=1)
    	      comp_res_all<-rbind(comp_res_all,comp_res)
    	      next
    	    } else {
    	      group1_r<-cor(scale(as.numeric(indep_pred[[group1]]$obs)), scale(indep_pred[[group1]]$pred))[1]
    	      group2_r<-cor(scale(as.numeric(indep_pred[[group2]]$obs)), scale(indep_pred[[group2]]$pred))[1]

    	      r_diff<-group1_r-group2_r

    	      group1_group2_r<-cor(scale(indep_pred[[group1]]$pred), scale(indep_pred[[group2]]$pred))

    	      r_diff_p<-paired.r(xy=group1_r, xz=group2_r, yz=group1_group2_r, n=length(scale(indep_pred[[group1]]$pred)), twotailed=T)$p[1]

    	      comp_res<-data.frame(Model_1=group1,
    	                           Model_2=group2,
    	                           Model_1_R=group1_r,
    	                           Model_2_R=group2_r,
    	                           R_diff=r_diff,
    	                           R_diff_pval=r_diff_p)
    	    }

    	    comp_res_all<-rbind(comp_res_all,comp_res)
    	 }
  	}

    # Write out the results
    write.table(comp_res_all, paste0(opt$output,'.pred_comp.txt'), col.names=T, row.names=F, quote=F)

    log_add(log_file = log_file, message = paste0('Model evaluation results saved as ',opt$output,'.pred_comp.txt.'))
  }

}

end.time <- Sys.time()
time.taken <- end.time - start.time
log_add(log_file = log_file, message = paste0('Analysis finished at ',as.character(end.time)))
log_add(log_file = log_file, message = paste0('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units')))
