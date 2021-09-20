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
		help="Prevalence of outcome in the general population [optional]"),
make_option("--out", action="store", default=NA, type='character',
		help="Prefix for output files [required]"),
make_option("--save_group_model", action="store", default=F, type='logical',
		help="Save group models for external validation [optional]"),
make_option("--assoc", action="store", default=T, type='logical',
		help="Perform association analysis between each predictor and outcome [optional]"),
make_option("--n_perm", action="store", default=1000, type='numeric',
		help="Number of permutations for model comparison [optional]"),
make_option("--compare_predictors", action="store", default=F, type='logical',
    help="Option to assign each predictor to own group [optional]"),
make_option("--eval_only", action="store", default=F, type='logical',
    help="Option to evaluate each model without comparison [optional]"),
make_option("--test", action="store", default=F, type='logical',
    help="Option to subset data for a test run [optional]"),
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
suppressMessages(library(psych))
suppressMessages(library(MASS))
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
if(opt$compare_predictors == F){
	if(!is.null(predictors_list$group)){
		if(length(unique(predictors_list$group)) != 1){
		  opt$group_info<-T
			opt$model_comp<-T
			predictors_list$group<-gsub("[[:punct:]]", ".",predictors_list$group)
			sink(file = paste(opt$out,'.log',sep=''), append = T)
			cat('Predictors file contains group information so model comparisons will be performed.\n')
			sink()
		} else {
		  opt$group_info<-T
			opt$model_comp<-F
			sink(file = paste(opt$out,'.log',sep=''), append = T)
			cat('Predictors file does not contain group information so model comparisons will not be performed.\n')
			sink()
		}
	} else {
	  opt$group_info<-F
		opt$model_comp<-F
		sink(file = paste(opt$out,'.log',sep=''), append = T)
		cat('Predictors file does not contain group information so model comparisons will not be performed.\n')
		sink()
	}
} else {
    if(!is.null(predictors_list$group)){
      opt$group_info<-T
      if(length(unique(predictors_list$group)) != 1){
        predictors_list$group<-gsub("[[:punct:]]", ".",predictors_list$group)
        sink(file = paste(opt$out,'.log',sep=''), append = T)
        cat('Predictors file contains group information so model comparisons will be performed.\n')
        sink()
      } else {
        sink(file = paste(opt$out,'.log',sep=''), append = T)
        cat('Predictors file does not contain group information.\n')
        sink()
      }
    } else {
      opt$group_info<-F
      sink(file = paste(opt$out,'.log',sep=''), append = T)
      cat('Predictors file does not contain group information.\n')
      sink()
    }
  
    sink(file = paste(opt$out,'.log',sep=''), append = T)
    cat('Each predictor will be assigned to a seperate group and model comparisons will be performed.\n')
    sink()
    
		opt$model_comp<-T
}

###########
# Read in the phenotypic data
###########

Outcome<-fread(opt$pheno)
Outcome<-Outcome[complete.cases(Outcome),]
names(Outcome)<-c('FID','IID','Outcome_var')

if(opt$test == T){
  Outcome<-Outcome[sample(1:dim(Outcome)[1],1000),]
}

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
	
	rm(keep_file)
	gc()
	
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

	if(opt$compare_predictors == T){
		# Create the object predictors_list
		predictors_list<-data.frame(predictor=names(Predictors)[-1], group=names(Predictors)[-1])
		predictors_list$group<-gsub("[[:punct:]]", ".",predictors_list$group)
		# Add group name to each predictor
		names(Predictors)[-1]<-paste0('Group_',names(Predictors)[-1],'.',names(Predictors)[-1])
	}
	
	sink(file = paste(opt$out,'.log',sep=''), append = T)
	cat('Predictors file contains',dim(Predictors)[2]-1,'predictors with sufficient data.\n')
	cat('Predictors file contains',dim(Predictors)[1],'individuals with complete data for remaining predictors.\n')
	sink()

}

predictors_list_new<-predictors_list

if(opt$compare_predictors == T){
  if(opt$group_info == T){
    # Create the object predictors_list
    predictors_list_new<-rbind(predictors_list_new,data.frame(predictor=names(Predictors)[-1], group=names(Predictors)[-1]))
    predictors_list_new$group<-gsub("[[:punct:]]", ".",predictors_list_new$group)
    predictors_list_new$group<-gsub('Group.','',predictors_list_new$group)
  }
}

if(opt$test == T & dim(predictors_list_new)[1] > 10){
  predictors_list_new<-predictors_list_new[1:10,]
  Predictors<-Predictors[,1:11]
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
				                            Estimate=NA,
				                            SE=NA,
				                            OR=NA,
				                            LowCI=NA,
				                            HighCI=NA,
				                            P=NA,
				                            AUC=NA,
				                            N=length(Outcome_Predictors_y),
				                            Ncas=sum(Outcome_Predictors_y == 'CASE'),
				                            Ncon=sum(Outcome_Predictors_y == 'CONTROL'),
				                            Obs_R2=NA)
				
				Assoc_res_temp
			} else {
			  mod<-glm(Outcome_Predictors_y ~ scale(Outcome_Predictors_x[,i]), family=opt$family)
			  prob<-predict(mod)
			  obs_r2<-cor(prob, as.numeric(Outcome_Predictors_y))^2
			  sum_mod<-summary(mod)
				roc_auc<-roc(Outcome_Predictors_y ~ prob)$auc
				cis<-exp(confint.default(mod))
				Assoc_res_temp<-data.frame(	Predictor=names(Outcome_Predictors_x)[i],
			                              Estimate=coef(sum_mod)[2,1],
			                              SE=coef(sum_mod)[2,2],
			                              OR=exp(coef(sum_mod)[2,1]),
			                              LowCI=cis[2,1],
			                              HighCI=cis[2,2],
				                            P=coef(sum_mod)[2,4],
				                            AUC=roc_auc,
				                            N=length(Outcome_Predictors_y),
				                            Ncas=sum(Outcome_Predictors_y == 'CASE'),
				                            Ncon=sum(Outcome_Predictors_y == 'CONTROL'),
				                            Obs_R2=obs_r2)
			  Assoc_res_temp
			}
		}
		# Convert R2 to liability scale
		Assoc_res$Liab_R2<-h2l_R2(opt$outcome_pop_prev, Assoc_res$Obs_R2, sum(Outcome_Predictors_y == 'CASE')/length(Outcome_Predictors_y))
	} else {
		Assoc_res<-foreach(i=1:dim(Outcome_Predictors_x)[2], .combine=rbind) %dopar% {
			if(var(Outcome_Predictors_x[,i]) == 0){
				Assoc_res_temp<-data.frame(	Predictor=names(Outcome_Predictors_x)[i],
																		BETA=NA,
																		SE=NA,
																		P=NA,
																		N=length(Outcome_Predictors_y),
																		Obs_R2=NA)
				
				Assoc_res_temp
			} else {
				mod<-glm(scale(Outcome_Predictors_y) ~ scale(Outcome_Predictors_x[,i]), family=opt$family)
			  sum_mod<-summary(mod)
			  Assoc_res_temp<-data.frame(	Predictor=names(Outcome_Predictors_x)[i],
											              BETA=coef(sum_mod)[2,1],
											              SE=coef(sum_mod)[2,2],
											              P=coef(sum_mod)[2,4],
											              N=length(Outcome_Predictors_y),
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


#print('451:')
#print(predictors_list_new)

# Create a model for each group, and using all predictors
predictors_list_new<-rbind(predictors_list_new,data.frame(predictor='All',
                                                          group='.'))
models<-list()
indep_pred<-list()

if(opt$model_comp == T){
	# Build glmnet using each group of predictors at a time
	for(group in unique(predictors_list_new$group)){
		# Subset predictor in the group
		print(group)
		Outcome_Predictors_train_x_group<-Outcome_Predictors_train_x[grepl(paste0('Group_',group,'\\.','|','Group_',group,'$'), names(Outcome_Predictors_train_x))]
		Outcome_Predictors_test_x_group<-Outcome_Predictors_test_x[grepl(paste0('Group_',group,'\\.','|','Group_',group,'$'), names(Outcome_Predictors_test_x))]
		
		# Rename group to 'all' if == '.'
		if(group == '.'){
		  group<-'All'
		  Outcome_Predictors_train_x_group<-Outcome_Predictors_train_x
		  Outcome_Predictors_test_x_group<-Outcome_Predictors_test_x
		}
		
		# If there is only one predictor, use glm
		if(dim(Outcome_Predictors_train_x_group)[2] > 1){
				model<- train(y=Outcome_Predictors_train_y, x=Outcome_Predictors_train_x_group, trControl=trainControl(method="cv", seeds=seeds, number=opt$n_fold, classProbs=T, savePredictions = 'final'), method="glmnet", family=opt$family)
		} else {
				model<- train(y=Outcome_Predictors_train_y, x=cbind(0,Outcome_Predictors_train_x_group), trControl=trainControl(method="cv", seeds=seeds, number=opt$n_fold, classProbs=T, savePredictions = 'final'), method="glm", family=opt$family)
		}
		
		if(opt$family=='binomial'){
		  Cross_mod<-summary(lm(scale(as.numeric(model$pred$obs)) ~ scale(model$pred$CASE)))
		  Cross_log_mod<-glm(model$pred$obs ~ scale(model$pred$CASE),family=opt$family)
		  Cross_log<-summary(Cross_log_mod)
		  Cross_cis<-exp(confint.default(Cross_log_mod))
		  Cross_auc<-roc(model$pred$obs ~ model$pred$CASE)$auc
		  Cross_LiabR2<-h2l_R2(opt$outcome_pop_prev, coef(Cross_mod)[2,1]^2, sum(Outcome_Predictors_train_y == 'CASE')/length(Outcome_Predictors_train_y))

			if(dim(Outcome_Predictors_train_x_group)[2] > 1){
				Indep_Pred<-predict(object = model$finalModel, newx = data.matrix(Outcome_Predictors_test_x_group), type = "response", s = model$finalModel$lambdaOpt)
			} else {
					tmp<-data.frame(cbind(0,Outcome_Predictors_test_x_group))
					names(tmp)[1]<-'0'
					Indep_Pred<-predict(object = model$finalModel, newdata = tmp, type = "response")
					rm(tmp)
			}
		  Indep_mod<-summary(lm(scale(as.numeric(Outcome_Predictors_test_y)) ~ scale(as.numeric(Indep_Pred))))
		  Indep_log_mod<-glm(Outcome_Predictors_test_y ~ scale(as.numeric(Indep_Pred)),family=opt$family)
		  Indep_log<-summary(Indep_log_mod)
		  Indep_cis<-exp(confint.default(Indep_log_mod))
		  Indep_auc<-roc(Outcome_Predictors_test_y ~ Indep_Pred)$auc
		  Indep_LiabR2<-h2l_R2(opt$outcome_pop_prev, coef(Indep_mod)[2,1]^2, sum(Outcome_Predictors_test_y == 'CASE')/length(Outcome_Predictors_test_y))

			Prediction_summary<-data.frame(	Model=paste0(group,'_group'),
			                                CrossVal_R=coef(Cross_mod)[2,1],
			                                CrossVal_R_SE=coef(Cross_mod)[2,2],
			                                CrossVal_OR=exp(coef(Cross_log)[2,1]),
			                                CrossVal_LowCI=Cross_cis[2,1],
			                                CrossVal_HighCI=Cross_cis[2,2],
			                                Cross_LiabR2=Cross_LiabR2,
			                                Cross_AUC=Cross_auc,
			                                CrossVal_pval=coef(Cross_mod)[2,4],
			                                CrossVal_N=length(model$pred$obs),
			                                CrossVal_Ncas=sum(model$pred$obs == 'CASE'),
			                                CrossVal_Ncon=sum(model$pred$obs == 'CONTROL'),
																			IndepVal_R=coef(Indep_mod)[2,1],
																			IndepVal_R_SE=coef(Indep_mod)[2,2],
																			IndepVal_OR=exp(coef(Indep_log)[2,1]),
																			IndepVal_LowCI=Indep_cis[2,1],
																			IndepVal_HighCI=Indep_cis[2,2],
																			Indep_LiabR2=Indep_LiabR2,
																			Indep_AUC=Indep_auc,
																			IndepVal_pval=coef(Indep_mod)[2,4],
																			IndepVal_N=length(Outcome_Predictors_test_y),
																			IndepVal_Ncas=sum(Outcome_Predictors_test_y == 'CASE'),
																			IndepVal_Ncon=sum(Outcome_Predictors_test_y == 'CONTROL'))
			
		} else {
			Cross_mod<-summary(lm(scale(model$pred$obs) ~ scale(model$pred$pred)))
			if(dim(Outcome_Predictors_train_x_group)[2] > 1){
				Indep_Pred<-predict(object = model$finalModel, newx = data.matrix(Outcome_Predictors_test_x_group), type = "response", s = model$finalModel$lambdaOpt)
			} else {
					tmp<-data.frame(cbind(0,Outcome_Predictors_test_x_group))
					names(tmp)[1]<-'0'
					Indep_Pred<-predict(object = model$finalModel, newdata = tmp, type = "response")
					rm(tmp)
			}
			Indep_mod<-summary(lm(scale(Outcome_Predictors_test_y) ~ scale(Indep_Pred)))

			Prediction_summary<-data.frame(	Model=paste0(group,'_group'),
																			CrossVal_R=coef(Cross_mod)[2,1],
																			CrossVal_R_SE=coef(Cross_mod)[2,2],
																			CrossVal_pval=coef(Cross_mod)[2,4],
																			CrossVal_N=length(model$pred$obs),
																			IndepVal_R=coef(Indep_mod)[2,1],
																			IndepVal_R_SE=coef(Indep_mod)[2,2],
																			IndepVal_pval=coef(Indep_mod)[2,4],
																			IndepVal_N=length(Outcome_Predictors_test_y))
																			
		}

		rm(Outcome_Predictors_train_x_group,Outcome_Predictors_test_x_group)
		gc()
		
        
        #print('557')
        #print(Prediction_summary_all)
        #print(colnames(Prediction_summary_all))
        #print(colnames(Prediction_summary))
		Prediction_summary_all<-rbind(Prediction_summary_all,Prediction_summary)
	
		if(opt$save_group_model == T){
			################
			# Save the model for use in external samples
			################                                

			saveRDS(model$finalModel, paste0(opt$out,'.',group,'_group.model.rds'))

			sink(file = paste(opt$out,'.log',sep=''), append = T)
			cat(group,' model saved as ',opt$out,'.',group,'_group.model.rds.\n',sep='')
			sink()
		}
		
		if(opt$eval_only == F){
		  # Rename model object for comparison between models
		  models[[group]]<-model
		  indep_pred[[group]]<-Indep_Pred
		} else {
		  rm(model,Indep_Pred)
		  gc()
		}
		
	}
  # Write out the results
  write.table(Prediction_summary_all, paste0(opt$out,'.pred_eval.txt'), col.names=T, row.names=F, quote=F)
  
  sink(file = paste(opt$out,'.log',sep=''), append = T)
  cat('Model evaluation results saved as ',opt$out,'.pred_eval.txt.\n',sep='')
  sink()
}

if(opt$model_comp == T & opt$eval_only == F){
	###################
	# Compare predictive utiliy of the different models
	###################
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
  	      if(opt$family=='binomial'){
      	    group1_r_Cross<-cor(as.numeric(models[[group1]]$pred$obs),models[[group1]]$pred$CASE)[1]
      			group2_r_Cross<-cor(as.numeric(models[[group2]]$pred$obs),models[[group2]]$pred$CASE)[1]
  
      			group1_r_Indep<-cor(as.numeric(Outcome_Predictors_test_y),indep_pred[[group1]])[1]
    			  group2_r_Indep<-cor(as.numeric(Outcome_Predictors_test_y),indep_pred[[group2]])[1]
  	      } else {
  	        group1_r_Cross<-cor(models[[group1]]$pred$obs,models[[group1]]$pred$pred)[1]
  	        group2_r_Cross<-cor(models[[group2]]$pred$obs,models[[group2]]$pred$pred)[1]

  	        group1_r_Indep<-cor(Outcome_Predictors_test_y,indep_pred[[group1]])[1]
  	        group2_r_Indep<-cor(Outcome_Predictors_test_y,indep_pred[[group2]])[1]
  	      }
  	      
  	      comp_res<-data.frame(Model_1=group1,
  	                           Model_2=group2,
  	                           Model_1_Cross_R=group1_r_Cross,
  	                           Model_2_Cross_R=group2_r_Cross,
  	                           Cross_R_diff=0,
  	                           Cross_R_diff_pval=1,
  	                           Model_1_Indep_R=group1_r_Indep,
  	                           Model_2_Indep_R=group2_r_Indep,
  	                           Indep_R_diff=0,
  	                           Indep_R_diff_pval=1)
  	      comp_res_all<-rbind(comp_res_all,comp_res)
  	      next
  	    }
  
  	  if(opt$family=='binomial'){
  			group1_r_Cross<-cor(as.numeric(models[[group1]]$pred$obs),models[[group1]]$pred$CASE)[1]
  			group2_r_Cross<-cor(as.numeric(models[[group2]]$pred$obs),models[[group2]]$pred$CASE)[1]
  			r_Cross_diff<-group1_r_Cross-group2_r_Cross
  			
  			group1_group2_r_Cross<-cor(models[[group1]]$pred$CASE, models[[group2]]$pred$CASE)
  			
  			r_Cross_diff_p<-paired.r(xy=group1_r_Cross, xz=group2_r_Cross, yz=group1_group2_r_Cross, n=length(models[[group1]]$pred$CASE), twotailed=T)$p[1]
  			
  			group1_r_Indep<-cor(as.numeric(Outcome_Predictors_test_y),indep_pred[[group1]])[1]
  			group2_r_Indep<-cor(as.numeric(Outcome_Predictors_test_y),indep_pred[[group2]])[1]
  			r_Indep_diff<-group1_r_Indep-group2_r_Indep
  			
  			group1_group2_r_Indep<-cor(indep_pred[[group1]], indep_pred[[group2]])
  			
   			r_Indep_diff_p<-paired.r(xy=group1_r_Indep, xz=group2_r_Indep, yz=group1_group2_r_Indep, n=length(indep_pred[[group2]]), twotailed=T)$p[1]
  	  } else {
  	    group1_r_Cross<-cor(models[[group1]]$pred$obs,models[[group1]]$pred$pred)[1]
  	    group2_r_Cross<-cor(models[[group2]]$pred$obs,models[[group2]]$pred$pred)[1]
  	    r_Cross_diff<-group1_r_Cross-group2_r_Cross
  	    
  	    group1_group2_r_Cross<-cor(models[[group1]]$pred$pred, models[[group2]]$pred$pred)
  	    
  	    r_Cross_diff_p<-paired.r(xy=group1_r_Cross, xz=group2_r_Cross, yz=group1_group2_r_Cross, n=length(models[[group1]]$pred$pred), twotailed=T)$p[1]
  	    
  	    group1_r_Indep<-cor(Outcome_Predictors_test_y,indep_pred[[group1]])[1]
  	    group2_r_Indep<-cor(Outcome_Predictors_test_y,indep_pred[[group2]])[1]
  	    r_Indep_diff<-group1_r_Indep-group2_r_Indep
  	    
  	    group1_group2_r_Indep<-cor(indep_pred[[group1]], indep_pred[[group2]])
  	    
  	    r_Indep_diff_p<-paired.r(xy=group1_r_Indep, xz=group2_r_Indep, yz=group1_group2_r_Indep, n=length(indep_pred[[group2]]), twotailed=T)$p[1]
  	  }

  			comp_res<-data.frame(Model_1=group1,
  														Model_2=group2,
  														Model_1_Cross_R=group1_r_Cross,
  														Model_2_Cross_R=group2_r_Cross,
  														Cross_R_diff=r_Cross_diff,
  														Cross_R_diff_pval=r_Cross_diff_p,
  														Model_1_Indep_R=group1_r_Indep,
  														Model_2_Indep_R=group2_r_Indep,
  														Indep_R_diff=r_Indep_diff,
  														Indep_R_diff_pval=r_Indep_diff_p)
  	    comp_res_all<-rbind(comp_res_all,comp_res)		
  	 }
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

