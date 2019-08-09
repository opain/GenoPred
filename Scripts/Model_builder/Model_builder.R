#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--pheno", action="store", default=NA, type='character',
		help="File containing phenotypic data [required]"),
make_option("--predictors", action="store", default=NA, type='character',
		help="File listing files containing predictors, with a groups column for model comparison [required]"),
make_option("--n_fold", action="store", default=10, type='numeric',
		help="Number of folds in cross validation [optional]"),
make_option("--n_core", action="store", default=10, type='numeric',
		help="Number of cores for parallel computing [optional]"),
make_option("--keep", action="store", default=NA, type='character',
		help="File containing list of individuals to include in analysis [optional]"),
make_option("--internal_validation_prop", action="store", default=0.2, type='numeric',
		help="Proportion of data that should be used for internal validation [optional]"),
make_option("--outcome_pop_prev", action="store", default=NA, type='numeric',
		help="Prevelance of outcome in the general population [optional]"),
make_option("--out", action="store", default=NA, type='character',
		help="Prefix for output files [required]"),
make_option("--save_model", action="store", default=T, type='logical',
		help="Save final model for extranl validation [optional]"),
make_option("--assoc", action="store", default=T, type='logical',
		help="Perform assocaition analysis between each predictor and outcome [optional]"),
make_option("--n_perm", action="store", default=1000, type='numeric',
		help="Number of permutations for model comparison [optional]"),
make_option("--pred_miss", action="store", default=0.1, type='numeric',
		help="Proportion of missing values allowed in predictor [optional]")		
)

opt = parse_args(OptionParser(option_list=option_list))

sink(file = paste(opt$out,'.log',sep=''), append = F)
cat(
'#################################################################
# Model_builder.R
# For questions contact Oliver Pain (oliver.pain@kcl.ac.uk)
#################################################################
Analysis started at',as.character(start.time),'
Options are:\n')

cat('Options are:\n')
print(opt)
cat('Analysis started at',as.character(start.time),'\n')
sink()

suppressMessages(library(data.table))
suppressMessages(library(glmnet))
suppressMessages(library(doMC))
suppressMessages(library(caret))
suppressMessages(library(pROC))
suppressMessages(library(verification))
registerDoMC(opt$n_core)

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

# Read in the predictors file
predictors_list<-data.frame(fread(opt$predictors))
names(predictors_list)[1]<-'predictor'

# Determine if there is a group column, and if there is more than 1 group.
if(!is.null(predictors_list$group)){
	if(length(unique(predictors_list$group)) != 1){
		opt$model_comp<-T
		predictors_list$group<-gsub('-','.',predictors_list$group)
		sink(file = paste(opt$out,'.log',sep=''), append = T)
		cat('Predictors file contains group information so model comparisons will be performed.\n')
		sink()
	} else {
		opt$model_comp<-F
		sink(file = paste(opt$out,'.log',sep=''), append = T)
		cat('Predictors file does not contain group information so model comparisons will not be performed.\n')
		sink()
	}
} else {
	opt$model_comp<-F
	sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat('Predictors file does not contain group information so model comparisons will not be performed.\n')
	sink()
}

###########
# Read in the phenotypic data
###########

Outcome<-fread(opt$pheno)
Outcome<-Outcome[complete.cases(Outcome),]
names(Outcome)[3]<-'Outcome_var'

# Create a column that combines FID and IID
Outcome$IID<-paste0(Outcome$FID,':',Outcome$IID)
Outcome$FID<-NULL

sink(file = paste(opt$out,'.log',sep=''), append = T)
cat('Phenotype file contains',dim(Outcome)[1],'individuals with complete data.\n')
sink()

# Determine whether outcome is binary or continuous and format accordingly
if(dim(unique(Outcome[,2]))[1] > 2){
	opt$family<-'gaussian'
}
if(dim(unique(Outcome[,2]))[1] == 2){
	opt$family<-'binomial'
	Outcome$Outcome_var<-factor(Outcome$Outcome_var, labels=c('CONTROL','CASE'))
}

sink(file = paste(opt$out,'.log',sep=''), append = T)
	if(opt$family == 'binomial'){
		cat('Phenotype is binary.\n')
	} else {
		cat('Phenotype is quantitative.\n')
	}
sink()

if(!is.na(opt$keep)){
	############
	# Extract individuals in the keep file
	############

	# Read in keep file
	keep_file<-fread(opt$keep)
	names(keep_file)[1:2]<-c('FID','IID')
	keep_file$IID<-paste0(keep_file$FID,':',keep_file$IID)
	keep_file$FID<-NULL
	# Extract keep indviduls from the phenotypic data
	Outcome<-Outcome[(Outcome$IID %in% keep_file$IID),]
	
	sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat('Phenotype file contains ',dim(Outcome)[1],' individuals after extraction of individuals in ',opt$keep,'.\n',sep='')
	sink()
} 
 
############
# Read in the predictor variables
############

predictor_merger<-function(x,y){
	return(merge(x,y,by='IID'))
}

if(dim(predictors_list)[1]>1){
		Predictors<-data.frame(foreach(k=1:dim(predictors_list)[1], .combine=predictor_merger) %dopar% {
			
			Predictors_temp<-fread(predictors_list$predictor[k])

			if(names(Predictors_temp)[1] == 'FID' & names(Predictors_temp)[2] == 'IID'){
				Predictors_temp$IID<-paste0(Predictors_temp$FID,':',Predictors_temp$IID)
				Predictors_temp$FID<-NULL
			} else {
				names(Predictors_temp)[1]<-'IID'
				Predictors_temp$IID<-paste0(Predictors_temp$IID,':',Predictors_temp$IID)
			}	
		
			# Remove variables with > opt$pred_miss missing values 
			Predictors_temp <- Predictors_temp[,colSums(is.na(Predictors_temp))/nrow(Predictors_temp) < opt$pred_miss, with=F]

			# Remove individuals with any missing data
			Predictors_temp<-Predictors_temp[complete.cases(Predictors_temp),]

			# Update column names to avoid duplciate column names between predictor files
			names(Predictors_temp)[-1]<-paste0('PredFile',k,'.',names(Predictors_temp)[-1])
			
			if(opt$model_comp == T){
					# Add group name to each predictor
					names(Predictors_temp)[-1]<-paste0('Group_',predictors_list$group[k],'.',names(Predictors_temp)[-1])
			}
			
			sink(file = paste(opt$out,'.log',sep=''), append = T)
			cat('Predictors file',k,'contains',dim(Predictors_temp)[2]-1,'predictors with sufficient data.\n')
			cat('Predictors file',k,'contains',dim(Predictors_temp)[1],'individuals with complete data for remaining predictors.\n')
			sink()
			
			Predictors_temp	
		})

sink(file = paste(opt$out,'.log',sep=''), append = T)
cat('After merging the',length(predictors_list),'Predictors files,', dim(Predictors)[2]-1,'predictors remain.\n')
cat('After merging the',length(predictors_list),'Predictors files,', dim(Predictors)[1],'individuals remain.\n')
sink()

} else {
		
	# Read in the predictor variables
	Predictors<-fread(predictors_list$predictor[1])

	if(names(Predictors)[1] == 'FID' & names(Predictors)[2] == 'IID'){
		Predictors$IID<-paste0(Predictors$FID,':',Predictors$IID)
		Predictors$FID<-NULL
	} else {
		names(Predictors)[1]<-'IID'
		Predictors$IID<-paste0(Predictors$IID,':',Predictors$IID)
	}

	# Remove variables with > opt$pred_miss missing values 
	Predictors <- Predictors[,colSums(is.na(Predictors))/nrow(Predictors) < opt$pred_miss, with=F]

	# Remove individuals with any missing data
	Predictors<-Predictors[complete.cases(Predictors),]

	sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat('Predictors file contains',dim(Predictors)[2]-1,'predictors with sufficient data.\n')
	cat('Predictors file contains',dim(Predictors)[1],'individuals with complete data for remaining predictors.\n')
	sink()

}
###########
# Merge the phenotype and predictor variables
###########
Outcome_Predictors <- data.frame(merge(Outcome,Predictors, by='IID'))

rm(Outcome,Predictors)

sink(file = paste(opt$out,'.log',sep=''), append = T)
cat(dim(Outcome_Predictors)[1],'individuals have both phenotypic and predictor data.\n')
sink()

# Report the size of the combined outcome and predictor data
sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat("Data to be carrried foward is ",format(object.size(Outcome_Predictors), units='auto'),".\n",sep='')
sink()

if(opt$assoc == T){

	############
	# Test association between Outcome and each variable in Predictors
	############

	sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat('Performing association analysis with each predictor...')
	sink()

	Outcome_Predictors_y<-Outcome_Predictors$Outcome_var
	Outcome_Predictors_x<-Outcome_Predictors[-1:-2]

	Assoc_res<-NULL

	if(opt$family == 'binomial'){
		Assoc_res<-foreach(i=1:dim(Outcome_Predictors_x)[2], .combine=rbind) %dopar% {
			if(var(Outcome_Predictors_x[,i]) == 0){
				Assoc_res_temp<-data.frame(	Predictor=names(Outcome_Predictors_x)[i],
																		BETA=NA,
																		SE=NA,
																		P=NA,
																		Obs_R2=NA)
				
				Assoc_res_temp
			} else {
			  mod<-glm(Outcome_Predictors_y ~ scale(Outcome_Predictors_x[,i]), family=opt$family)
				obs_r2<-cor(predict(mod), as.numeric(Outcome_Predictors_y))^2
			  sum_mod<-summary(mod)
				prob<-predict(mod,type=c("response"))
			  Assoc_res_temp<-data.frame(	Predictor=names(Outcome_Predictors_x)[i],
											              BETA=coef(sum_mod)[2,1],
											              SE=coef(sum_mod)[2,2],
											              P=coef(sum_mod)[2,4],
											              Obs_R2=obs_r2)
			  Assoc_res_temp
			}
		}
		# Convert Nagelkerke R2 to liability scale
		Assoc_res$Liab_R2<-h2l_R2(opt$outcome_pop_prev, Assoc_res$Obs_R2, sum(Outcome_Predictors_y == 'CASE')/length(Outcome_Predictors_y))
	} else {
		Assoc_res<-foreach(i=1:dim(Outcome_Predictors_x)[2], .combine=rbind) %dopar% {
			if(var(Outcome_Predictors_x[,i]) == 0){
				Assoc_res_temp<-data.frame(	Predictor=names(Outcome_Predictors_x)[i],
																		BETA=NA,
																		SE=NA,
																		P=NA,
																		Obs_R2=NA)
				
				Assoc_res_temp
			} else {
				mod<-glm(scale(Outcome_Predictors_y) ~ scale(Outcome_Predictors_x[,i]), family=opt$family)
			  sum_mod<-summary(mod)
			  Assoc_res_temp<-data.frame(	Predictor=names(Outcome_Predictors_x)[i],
											              BETA=coef(sum_mod)[2,1],
											              SE=coef(sum_mod)[2,2],
											              P=coef(sum_mod)[2,4],
											              Obs_R2=coef(sum_mod)[2,1]^2)
				Assoc_res_temp
			}
		}
	}
	sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat('Done!\n')
	sink()

	# Write out the results
	write.table(Assoc_res, paste0(opt$out,'.assoc.txt'), col.names=T, row.names=F, quote=F)

	sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat('Predictor association results saved as ',opt$out,'.assoc.txt.\n',sep='')
	sink()
}

############
# Build and evaluate models using predictors together
############

# Split the sample into training and test data (4/5 training, 1/5 test)
set.seed(1)
Outcome_Predictors_train_ind <- sample(seq_len(dim(Outcome_Predictors)[1]), size = floor((1-opt$internal_validation_prop) * dim(Outcome_Predictors)[1]))

Outcome_Predictors_train <- Outcome_Predictors[Outcome_Predictors_train_ind, ]
Outcome_Predictors_test <- Outcome_Predictors[-Outcome_Predictors_train_ind, ]

Outcome_Predictors_train_y<-Outcome_Predictors_train$Outcome_var
Outcome_Predictors_train_x<-Outcome_Predictors_train[-1:-2]

Outcome_Predictors_test_y<-Outcome_Predictors_test$Outcome_var
Outcome_Predictors_test_x<-Outcome_Predictors_test[-1:-2]

rm(Outcome_Predictors,Outcome_Predictors_train,Outcome_Predictors_test)

if(opt$family == 'binomial'){
	sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat("Training data contains ", length(Outcome_Predictors_train_y)," individuals (",sum(Outcome_Predictors_train_y == 'CASE')," cases and ",sum(Outcome_Predictors_train_y == 'CONTROL')," controls).\n",sep='')
	sink()
	sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat("Test data contains ", length(Outcome_Predictors_test_y)," individuals (",sum(Outcome_Predictors_test_y == 'CASE')," cases and ",sum(Outcome_Predictors_test_y == 'CONTROL')," controls).\n",sep='')
	sink()
} else {
	sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat("Training data contains ", length(Outcome_Predictors_train_y)," individuals.\n",sep='')
	sink()
	sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat("Test data contains ", length(Outcome_Predictors_test_y)," individuals.\n",sep='')
	sink()
}

# Report the size of the predictors data.frame
sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat("Outcome data in training sample is ",format(object.size(Outcome_Predictors_train_y), units='auto'),".\n",sep='')
	cat("Predictor data in training sample is ",format(object.size(Outcome_Predictors_train_x), units='auto'),".\n", sep='')
sink()

# Create a variable containing seeds for cross validation. (this must be set to reproducible results)
set.seed(1)
seeds <- vector(mode = "list", length = opt$n_fold+1)
for(i in 1:(opt$n_fold)) seeds[[i]]<- sample.int(n=1000, 10)
seeds[[opt$n_fold+1]]<-sample.int(n=1000, 1)

Prediction_summary_all<-NULL

if(opt$model_comp == T){
	# Build glmnet using each group of predictors at a time
	for(group in unique(predictors_list$group)){
		# Subset predictor in the group
		Outcome_Predictors_train_x_group<-Outcome_Predictors_train_x[grepl(paste0('Group_',group),names(Outcome_Predictors_train_x))]
		Outcome_Predictors_test_x_group<-Outcome_Predictors_test_x[grepl(paste0('Group_',group),names(Outcome_Predictors_test_x))]
		
		# If there is only one predictor, add empty variable so it runs
		if(dim(Outcome_Predictors_train_x_group)[2] > 1){
			model<- train(y=Outcome_Predictors_train_y, x=Outcome_Predictors_train_x_group, trControl=trainControl(method="cv", seeds=seeds, number=opt$n_fold, classProbs=T, savePredictions = 'final'), method="glmnet", family=opt$family)
		} else {
			model<- train(y=Outcome_Predictors_train_y, x=cbind(0,Outcome_Predictors_train_x_group), trControl=trainControl(method="cv", seeds=seeds, number=opt$n_fold, classProbs=T, savePredictions = 'final'), method="glmnet", family=opt$family)
		}
		
		if(opt$family=='binomial'){
			Cross_mod<-summary(lm(scale(as.numeric(model$pred$obs)) ~ scale(model$pred$CASE)))
			Cross_LiabR2<-h2l_R2(opt$outcome_pop_prev, coef(Cross_mod)[2,1]^2, sum(Outcome_Predictors_train_y == 'CASE')/length(Outcome_Predictors_train_y))

			if(dim(Outcome_Predictors_train_x_group)[2] > 1){
				Indep_Pred<-predict(object = model$finalModel, newx = data.matrix(Outcome_Predictors_test_x_group), type = "response", s = model$finalModel$lambdaOpt)
			} else {
				Indep_Pred<-predict(object = model$finalModel, newx = data.matrix(cbind(0,Outcome_Predictors_test_x_group)), type = "response", s = model$finalModel$lambdaOpt)
			}
			Indep_mod<-summary(lm(scale(as.numeric(Outcome_Predictors_test_y)) ~ scale(as.numeric(Indep_Pred))))
			Indep_LiabR2<-h2l_R2(opt$outcome_pop_prev, coef(Indep_mod)[2,1]^2, sum(Outcome_Predictors_test_y == 'CASE')/length(Outcome_Predictors_test_y))

			Prediction_summary<-data.frame(	Model=paste0(group,'_group'),
																			CrossVal_R=coef(Cross_mod)[2,1],
																			CrossVal_R_SE=coef(Cross_mod)[2,2],
																			Cross_LiabR2=Cross_LiabR2,
																			CrossVal_pval=coef(Cross_mod)[2,4],
																			IndepVal_R=coef(Indep_mod)[2,1],
																			IndepVal_R_SE=coef(Indep_mod)[2,2],
																			Indep_LiabR2=Indep_LiabR2,
																			IndepVal_pval=coef(Indep_mod)[2,4])
			
			# Rename model object for comparison between models
			assign(paste0(group,"_model"),model)
			assign(paste0(group,"_Indep_Pred"),Indep_Pred)
		} else {
			Cross_mod<-summary(lm(scale(model$pred$obs) ~ scale(model$pred$pred)))
			if(dim(Outcome_Predictors_train_x_group)[2] > 1){
				Indep_Pred<-predict(object = model$finalModel, newx = data.matrix(Outcome_Predictors_test_x_group), type = "response", s = model$finalModel$lambdaOpt)
			} else {
				Indep_Pred<-predict(object = model$finalModel, newx = data.matrix(cbind(0,Outcome_Predictors_test_x_group)), type = "response", s = model$finalModel$lambdaOpt)
			}
			Indep_mod<-summary(lm(scale(Outcome_Predictors_test_y) ~ scale(Indep_Pred)))

			Prediction_summary<-data.frame(	Model=paste0(group,'_group'),
																			CrossVal_R=coef(Cross_mod)[2,1],
																			CrossVal_R_SE=coef(Cross_mod)[2,2],
																			CrossVal_pval=coef(Cross_mod)[2,4],
																			IndepVal_R=coef(Indep_mod)[2,1],
																			IndepVal_R_SE=coef(Indep_mod)[2,2],
																			IndepVal_pval=coef(Indep_mod)[2,4])
																			
			# Rename model object for comparison between models
			assign(paste0(group,"_model"),model)
			assign(paste0(group,"_Indep_Pred"),Indep_Pred)
		}
		
		Prediction_summary_all<-rbind(Prediction_summary_all,Prediction_summary)
	
		if(opt$save_model == T){
			################
			# Save the model for use in external samples
			################                                

			saveRDS(model$finalModel, paste0(opt$out,'.',group,'_group.model.rds'))

			sink(file = paste(opt$out,'.log',sep=''), append = T)
			cat(group,' model saved as ',opt$out,'.',group,'_group.model.rds.\n',sep='')
			sink()
		}
	}
}

#########
# Build glmnet using all predictors at once
#########
if(dim(Outcome_Predictors_train_x)[2] > 1){
	model<- train(y=Outcome_Predictors_train_y, x=Outcome_Predictors_train_x, trControl=trainControl(method="cv", seeds=seeds, number=opt$n_fold, classProbs=T, savePredictions = 'final'), method="glmnet", family=opt$family)
} else {
	model<- train(y=Outcome_Predictors_train_y, x=cbind(0,Outcome_Predictors_train_x), trControl=trainControl(method="cv", seeds=seeds, number=opt$n_fold, classProbs=T, savePredictions = 'final'), method="glmnet", family=opt$family)
}

if(opt$family=='binomial'){
	Cross_mod<-summary(lm(scale(as.numeric(model$pred$obs)) ~ scale(model$pred$CASE)))
	Cross_LiabR2<-h2l_R2(opt$outcome_pop_prev, coef(Cross_mod)[2,1]^2, sum(Outcome_Predictors_train_y == 'CASE')/length(Outcome_Predictors_train_y))

	
	if(dim(Outcome_Predictors_train_x)[2] > 1){
		Pred<-predict(object = model$finalModel, newx = data.matrix(Outcome_Predictors_test_x), type = "response", s = model$finalModel$lambdaOpt)
	} else {
		Pred<-predict(object = model$finalModel, newx = data.matrix(cbind(0,Outcome_Predictors_test_x)), type = "response", s = model$finalModel$lambdaOpt)
	}
	
	Indep_mod<-summary(lm(scale(as.numeric(Outcome_Predictors_test_y)) ~ scale(as.numeric(Pred))))

	Indep_LiabR2<-h2l_R2(opt$outcome_pop_prev, coef(Indep_mod)[2,1]^2, sum(Outcome_Predictors_test_y == 'CASE')/length(Outcome_Predictors_test_y))

	Prediction_summary<-data.frame(	Model='Full_model',
																	CrossVal_R=coef(Cross_mod)[2,1],
																	CrossVal_R_SE=coef(Cross_mod)[2,2],
																	Cross_LiabR2=Cross_LiabR2,
																	CrossVal_pval=coef(Cross_mod)[2,4],
																	IndepVal_R=coef(Indep_mod)[2,1],
																	IndepVal_R_SE=coef(Indep_mod)[2,2],
																	Indep_LiabR2=Indep_LiabR2,
																	IndepVal_pval=coef(Indep_mod)[2,4])
} else {
	Cross_mod<-summary(lm(scale(as.numeric(model$pred$obs)) ~ scale(model$pred$pred)))
	if(dim(Outcome_Predictors_train_x)[2] > 1){
		Pred<-predict(object = model$finalModel, newx = data.matrix(Outcome_Predictors_test_x), type = "response", s = model$finalModel$lambdaOpt)
	} else {
		Pred<-predict(object = model$finalModel, newx = data.matrix(cbind(0,Outcome_Predictors_test_x)), type = "response", s = model$finalModel$lambdaOpt)
	}
	
	Indep_mod<-summary(lm(scale(as.numeric(Outcome_Predictors_test_y)) ~ scale(Pred)))
	Prediction_summary<-data.frame(	Model='Full_model',
																	CrossVal_R=coef(Cross_mod)[2,1],
																	CrossVal_R_SE=coef(Cross_mod)[2,2],
																	CrossVal_pval=coef(Cross_mod)[2,4],
																	IndepVal_R=coef(Indep_mod)[2,1],
																	IndepVal_R_SE=coef(Indep_mod)[2,2],
																	IndepVal_pval=coef(Indep_mod)[2,4])
}

Prediction_summary_all<-rbind(Prediction_summary_all,Prediction_summary)
	
# Write out the results
write.table(Prediction_summary_all, paste0(opt$out,'.pred_eval.txt'), col.names=T, row.names=F, quote=F)

sink(file = paste(opt$out,'.log',sep=''), append = T)
cat('Model evaluation results saved as ',opt$out,'.pred_eval.txt.\n',sep='')
sink()

if(opt$save_model == T){
	################
	# Save the model for use in external samples
	################                                

	saveRDS(model$finalModel, paste0(opt$out,'.final_model.rds'))

	sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat('Final model saved as ',opt$out,'.final_model.rds.\n',sep='')
	sink()
}

if(opt$model_comp == T){
	###################
	# Compare predictive utiliy of the different models
	###################
	comp_res_all<-NULL
	for(group in unique(predictors_list$group)){
		if(opt$family=='binomial'){
			full_r_Cross<-cor(as.numeric(model$pred$obs),model$pred$CASE)[1]
			nest_r_Cross<-cor(as.numeric(get(paste0(group,"_model"))$pred$obs),get(paste0(group,"_model"))$pred$CASE)[1]
			r_Cross_diff<-abs(full_r_Cross)-abs(nest_r_Cross)
			
			r_Cross_diff_perm_all<-foreach(perm=1:opt$n_perm, .combine=c) %dopar% {
				full_r_Cross_perm<-cor(as.numeric(model$pred$obs),sample(model$pred$CASE))[1]
				nest_r_Cross_perm<-cor(as.numeric(get(paste0(group,"_model"))$pred$obs),sample(get(paste0(group,"_model"))$pred$CASE))[1]
				r_Cross_diff_perm<-abs(full_r_Cross_perm)-abs(nest_r_Cross_perm)
				r_Cross_diff_perm
			}
			
			r_Cross_diff_p<-as.character(sum(r_Cross_diff_perm_all >= r_Cross_diff)/length(r_Cross_diff_perm_all))
			if(r_Cross_diff_p == '0'){
					r_Cross_diff_p<-paste0('<',1/opt$n_perm)
			}
			
			full_r_Indep<-cor(as.numeric(Outcome_Predictors_test_y),Pred)[1]
			nest_r_Indep<-cor(as.numeric(Outcome_Predictors_test_y),get(paste0(group,"_Indep_Pred")))[1]
			r_Indep_diff<-abs(full_r_Indep)-abs(nest_r_Indep)
			
			r_Indep_diff_perm_all<-foreach(perm=1:opt$n_perm, .combine=c) %dopar% {
				full_r_Indep_perm<-cor(as.numeric(Outcome_Predictors_test_y),sample(Pred))[1]
				nest_r_Indep_perm<-cor(as.numeric(Outcome_Predictors_test_y),sample(get(paste0(group,"_Indep_Pred"))))[1]
				r_Indep_diff_perm<-abs(full_r_Indep_perm)-abs(nest_r_Indep_perm)
				r_Indep_diff_perm
			}
			
			r_Indep_diff_p<-as.character(sum(r_Indep_diff_perm_all >= r_Indep_diff)/length(r_Indep_diff_perm_all))
			if(r_Indep_diff_p == '0'){
					r_Indep_diff_p<-paste0('<',1/opt$n_perm)
			}
			
			comp_res<-data.frame(	Model_1='Full_model',
														Model_2=paste0(group,'_only_model'),
														Model_1_Cross_R=full_r_Cross,
														Model_2_Cross_R=nest_r_Cross,
														Cross_R_diff=r_Cross_diff,
														Cross_R_diff_pval=r_Cross_diff_p,
														Model_1_Indep_R=full_r_Indep,
														Model_2_Indep_R=nest_r_Indep,
														Indep_R_diff=r_Indep_diff,
														Indep_R_diff_pval=r_Indep_diff_p)
			} else {			
			full_r_Cross<-cor(model$pred$obs,model$pred$pred)[1]
			nest_r_Cross<-cor(get(paste0(group,"_model"))$pred$obs,get(paste0(group,"_model"))$pred$pred)[1]
			r_Cross_diff<-abs(full_r_Cross)-abs(nest_r_Cross)
			
			r_Cross_diff_perm_all<-foreach(perm=1:opt$n_perm, .combine=c) %dopar% {
				full_r_Cross_perm<-cor(model$pred$obs,sample(model$pred$pred))[1]
				nest_r_Cross_perm<-cor(get(paste0(group,"_model"))$pred$obs,sample(get(paste0(group,"_model"))$pred$pred))[1]
				r_Cross_diff_perm<-abs(full_r_Cross_perm)-abs(nest_r_Cross_perm)
				r_Cross_diff_perm
			}
			
			r_Cross_diff_p<-as.character(sum(r_Cross_diff_perm_all >= r_Cross_diff)/length(r_Cross_diff_perm_all))
			if(r_Cross_diff_p == '0'){
					r_Cross_diff_p<-paste0('<',1/opt$n_perm)
			}
			
			full_r_Indep<-cor(Outcome_Predictors_test_y,Pred)[1]
			nest_r_Indep<-cor(Outcome_Predictors_test_y,get(paste0(group,"_Indep_Pred")))[1]
			r_Indep_diff<-abs(full_r_Indep)-abs(nest_r_Indep)
			
			r_Indep_diff_perm_all<-foreach(perm=1:opt$n_perm, .combine=c) %dopar% {
				full_r_Indep_perm<-cor(Outcome_Predictors_test_y,sample(Pred))[1]
				nest_r_Indep_perm<-cor(Outcome_Predictors_test_y,sample(get(paste0(group,"_Indep_Pred"))))[1]
				r_Indep_diff_perm<-abs(full_r_Indep_perm)-abs(nest_r_Indep_perm)
				r_Indep_diff_perm
			}
			
			r_Indep_diff_p<-as.character(sum(r_Indep_diff_perm_all >= r_Indep_diff)/length(r_Indep_diff_perm_all))
			if(r_Indep_diff_p == '0'){
					r_Indep_diff_p<-paste0('<',1/opt$n_perm)
			}
			
			comp_res<-data.frame(	Model_1='Full_model',
														Model_2=paste0(group,'_only_model'),
														Model_1_Cross_R=full_r_Cross,
														Model_2_Cross_R=nest_r_Cross,
														Cross_R_diff=r_Cross_diff,
														Cross_R_diff_pval=r_Cross_diff_p,
														Model_1_Indep_R=full_r_Indep,
														Model_2_Indep_R=nest_r_Indep,
														Indep_R_diff=r_Indep_diff,
														Indep_R_diff_pval=r_Indep_diff_p)
			}
	comp_res_all<-rbind(comp_res_all,comp_res)		
	}

# Write out the results
write.table(comp_res_all, paste0(opt$out,'.pred_comp.txt'), col.names=T, row.names=F, quote=F)

sink(file = paste(opt$out,'.log',sep=''), append = T)
cat('Model evaluation results saved as ',opt$out,'.pred_comp.txt.\n',sep='')
sink()
}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$out,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
