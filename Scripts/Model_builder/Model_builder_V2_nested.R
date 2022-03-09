#!/usr/bin/Rscript
# This script was written by Oliver Pain whilst at King's College London University.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--pheno", action="store", default=NA, type='character',
		help="File containing phenotypic data [required]"),
make_option("--predictors", action="store", default=NA, type='character',
		help="File listing files containing predictors, with a groups column for model comparison [required]"),
make_option("--n_outer_fold", action="store", default=5, type='numeric',
    help="Number of folds in for outer cross-validation [optional]"),
make_option("--n_inner_fold", action="store", default=10, type='numeric',
    help="Number of folds for inner cross-validation [optional]"),
make_option("--n_core", action="store", default=1, type='numeric',
		help="Number of cores for parallel computing [optional]"),
make_option("--keep", action="store", default=NA, type='character',
		help="File containing list of individuals to include in analysis [optional]"),
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
make_option("--interaction", action="store", default=F, type='logical',
    help="Option to include interaction terms between predictors in each group [optional]"),
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
			col_keep<-T
			for(i in 2:ncol(Predictors_temp)){
			  col_keep<-c(col_keep, sum(!is.finite(Predictors_temp[[names(Predictors_temp)[i]]]) | is.na(Predictors_temp[[names(Predictors_temp)[i]]]))/nrow(Predictors_temp) < opt$pred_miss)
			}
			
			Predictors_temp <- Predictors_temp[,col_keep, with=F]

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
	col_keep<-T
	for(i in 2:ncol(Predictors_temp)){
	  col_keep<-c(col_keep, sum(!is.finite(Predictors_temp[[names(Predictors_temp)[i]]]) | is.na(Predictors_temp[[names(Predictors_temp)[i]]]))/nrow(Predictors_temp) < opt$pred_miss)
	}
	
	Predictors_temp <- Predictors_temp[,col_keep, with=F]
	
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

if(opt$model_comp == T){
############
# Build and evaluate models using predictors together
############

# In this version of Model_builder, we will use nested cross validation.
# We split the data into test and training samples, derive the model using 10-fold cross validation in the training sample, evaluate the model in the test sample, and then repeat this process until all parts of the sample have been used as the test sample.
# i.e. What I was doing before but repeated several time.

# Split the sample into n=opt$internal_validation_prop equal parts
set.seed(1)
nr<-dim(Outcome_Predictors)[1]
d<-sample(1:nr)

train.ext=createFolds(d,k=opt$n_outer_fold,returnTrain=TRUE)
test.ext=lapply(train.ext,function(x) (1:nr)[-x])

# Create a variable containing seeds for internal cross validation
set.seed(1)
seeds <- vector(mode = "list", length = opt$n_inner_fold+1)
for(i in 1:(opt$n_inner_fold)) seeds[[i]]<- sample.int(n=1000, 10)
seeds[[opt$n_inner_fold+1]]<-sample.int(n=1000, 1)

Prediction_summary_all<-NULL

# Create a model for each group, and using all predictors
predictors_list_new<-rbind(predictors_list_new,data.frame(predictor='All',
                                                          group='.'))

sink(file = paste(opt$out,'.log',sep=''), append = T)
cat("Initiating nested cross-validation...\n",sep='')
sink()

indep_pred<-list()
Prediction_summary_all<-NULL

if(opt$model_comp == T){
  # Build glmnet using each group of predictors at a time
  for(group in unique(predictors_list_new$group)){
    for(outer_val in 1:opt$n_outer_fold){
      print(outer_val)
        
      Outcome_Predictors_train_ind<-d[train.ext[[outer_val]]]
          
      Outcome_Predictors_train <- Outcome_Predictors[Outcome_Predictors_train_ind, ]
      Outcome_Predictors_test <- Outcome_Predictors[-Outcome_Predictors_train_ind, ]
      
      Outcome_Predictors_train_y<-Outcome_Predictors_train$Outcome_var
      Outcome_Predictors_train_x<-Outcome_Predictors_train[-1:-2]
      
      Outcome_Predictors_test_y<-Outcome_Predictors_test$Outcome_var
      Outcome_Predictors_test_x<-Outcome_Predictors_test[-1:-2]
      
  		# Subset predictor in the group
  		print(group)
  		Outcome_Predictors_train_x_group<-Outcome_Predictors_train_x[grepl(paste0('Group_',group,'\\.','|','Group_',group,'$'), names(Outcome_Predictors_train_x))]
  		Outcome_Predictors_test_x_group<-Outcome_Predictors_test_x[grepl(paste0('Group_',group,'\\.','|','Group_',group,'$'), names(Outcome_Predictors_test_x))]
		
  		# Rename group to 'all' if == '.'
  		if(group == '.'){
  		  group_name<-'All'
  		  Outcome_Predictors_train_x_group<-Outcome_Predictors_train_x
  		  Outcome_Predictors_test_x_group<-Outcome_Predictors_test_x
  		} else {
  		  group_name<-group
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
  write.table(Prediction_summary_all, paste0(opt$out,'.pred_eval.txt'), col.names=T, row.names=F, quote=F)
  
  sink(file = paste(opt$out,'.log',sep=''), append = T)
  cat('Model evaluation results saved as ',opt$out,'.pred_eval.txt.\n',sep='')
  sink()

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
  write.table(comp_res_all, paste0(opt$out,'.pred_comp.txt'), col.names=T, row.names=F, quote=F)
  
  sink(file = paste(opt$out,'.log',sep=''), append = T)
  cat('Model evaluation results saved as ',opt$out,'.pred_comp.txt.\n',sep='')
  sink()
}

}

end.time <- Sys.time()
time.taken <- end.time - start.time
sink(file = paste(opt$out,'.log',sep=''), append = T)
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
sink()
