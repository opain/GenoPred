# Model_builder_V2.R

This is script is used to evaluate and compare predictors, using an elastic net to model multiple predictors. Model_builder_V2.R uses 10-fold cross validation to tune elastic net hyperparameters and then tests the best model in an independent test set, avoiding overfitting due to hyperparameter optimisation. The predictive utility of different models are compared using the paired.r function from the psych package, which accounts for the correlation between models.

## Pre-requisites
The following software is required for the prediction pipeline:

* R packages:
```R
install.packages(c('data.table','glmnet','doMC','caret','pROC','verification','psych','MASS'))
```

## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --pheno | File containing phenotypic data [required] | NA |
| --predictors | File listing files containing predictors, with a group column for model comparison [required] | NA |
| --n_fold | Number of folds in cross validation [optional] | 10 |
| --n_core | Number of cores for parallel computing [optional] | 10 |
| --keep | File containing list of individuals to include in analysis [optional] | NA |
| --internal_validation_prop | Proportion of data that should be used for internal validation [optional] | 0.2 |
| --outcome_pop_prev | Prevalence of outcome in the general population [optional] | NA |
| --out | Prefix for output files [required] | NA |
| --save_group_model | Save group models for external validation [optional] | FALSE |
| --assoc | Perform association analysis between each predictor and outcome [optional] | TRUE |
| --n_perm | Number of permutations for model comparison [optional] | 1000 |
| --compare_predictors | Option to assign each predictor to own group [optional] | FALSE |
| --eval_only | Option to evaluate each model without comparison [optional] | FALSE |
| --test | Option to subset data for a test run [optional] | FALSE |
| --pred_miss | Proportion of missing values allowed in predictor [optional] | 0.1 |

## Input files

* --predictors: Should be a file with two columns labelled predictors and group. The predictors column specifies files containing predictor data, which should have columns FID, IID, and then on or more predictor variables. The group column specifies which predictors should be combined into a single model.
* --pheno: File containing outcome data with columns FID, IID, and then one outcome variable. If outcome is binary, logistic regression will be used. Otherwise, linear regression will be used.
* 

## Output files

The script outputs three for files:

* .assoc.txt: Univariate linear (continuous outcome)/logistic (binary outcome) regression results for each predictor
  * Linear model:
    * Predictor: Name of model
    * BETA: Effect size in units of standard deviation
    * SE: Standard error of BETA
    * P: P-value for association
    * Obs_R2: R-squared/variance explained on the observed scale
  * Logistic:
    * Predictor: Name of model
    * Estimate: Effect size with predictor in units of standard deviation
    * SE: Standard error of Estimate
    * OR: Odds Ratio
    * LowCI: 2.5% confidence interval for OR
    * HighCI: 97.5% confidence interval for OR
    * P: P-value for association
    * AUC: Area-under-the-ROC-curve
    * N: Sample size
    * Ncas: Number of cases
    * Ncon: Number of controls
    * Obs_R2: R-squared/variance explained on the observed scale
    * Liab_R2: R-squared/variance explained on the liability scale
* .pred_eval.txt: Predictive utility of each model
  * Linear:
    * Model: Model name
    * CrossVal_R: Pearson correlation between model predictions and outcome, as estimated using 10-fold cross validation
    * CrossVal_R_SE: Standard error of CrossVal_R
    * CrossVal_pval: p-value of CrossVal_R
    * CrossVal_N: Sample size in 10-fold cross validation
    * IndepVal_R: Pearson correlation between model predictions and outcome, as estimated using independent test set validation
    * IndepVal_R_SE: Standard error of IndepVal_R
    * IndepVal_pval: p-value of IndepVal_R
    * IndepVal_N: Sample size in independent test set validation
  * Logistic:
    * Model: Model name
    * CrossVal_R: Pearson correlation between model predictions and outcome, as estimated using 10-fold cross validation
    * CrossVal_R_SE: Standard error of CrossVal_R
    * CrossVal_OR: Odds ratio, as estimated using 10-fold cross validation
    * CrossVal_LowCI: 2.5% confidence interval of CrossVal_OR
    * CrossVal_HighCI: 97.5% confidence interval of CrossVal_OR
    * Cross_LiabR2: R-squared/variance explained on the liability scale, as estimated using 10-fold cross validation
    * Cross_AUC: Area-under-the-ROC-curve, as estimated using 10-fold cross validation
    * CrossVal_pval: p-value of CrossVal_R
    * CrossVal_N: Sample size in 10-fold cross validation
    * CrossVal_Ncas: Number of cases in 10-fold cross validation
    * CrossVal_Ncon: Number of controls in 10-fold cross validation
    * IndepVal_R: Pearson correlation between model predictions and outcome, as estimated using independent test set validation
    * IndepVal_R_SE: Standard error of IndepVal_R
    * IndepVal_OR: Odds ratio, as estimated using independent test set validation
    * IndepVal_LowCI: 2.5% confidence interval of IndepVal_OR
    * IndepVal_HighCI: 97.5% confidence interval of IndepVal_OR
    * Cross_LiabR2: R-squared/variance explained on the liability scale, as estimated using independent test set validation
    * Cross_AUC: Area-under-the-ROC-curve, as estimated using independent test set validation
    * IndepVal_pval: p-value of IndepVal_R
    * IndepVal_N: Sample size in independent test set validation
    * IndepVal_Ncas: Number of cases in independent test set validation
    * IndepVal_Ncon: Number of controls in independent test set validation
* .pred_comp.txt: Compares predictive utility of each model
  * Linear/Logistic:
    * Model_1: Name of model 1
    * Model_2: Name of model 2
    * Model_1_Cross_R: Pearson correlation between model 1 predictions and outcome, as estimated using 10-fold cross validation
    * Model_2_Cross_R: Pearson correlation between model 2 predictions and outcome, as estimated using 10-fold cross validation
    * Cross_R_diff: Model_1_Cross_R - Model_2_Cross_R
    * Cross_R_diff_pval: P-value of Cross_R_diff
    * Model_1_Indep_R: Pearson correlation between model 1 predictions and outcome, as estimated using independent test set validation
    * Model_2_Indep_R: Pearson correlation between model 2 predictions and outcome, as estimated using independent test set validation
    * Indep_R_diff: Model_1_Indep_R - Model_2_Indep_R
    * Indep_R_diff_pval: P-value of Indep_R_diff
* .log: Log file

## Examples
```sh
Rscript Model_builder/Model_builder_V2.R \
    --pheno phenotype_data.txt \
    --keep indiv_keep.txt \
    --out example \
    --n_core 2 \
    --compare_predictors F \
    --assoc T \
    --outcome_pop_prev 0.3 \
    --predictors predictors_txt
```

***

# Model_builder_V2_nested.R

Alsmost exactly the same as Model_builder_V2.R except it uses nested 10-fold cross validation, iteratively using all parts of the sample as the independent test set. Therefore, the output only contains one set of predictive utility results for each model, rather than one set for cross validation and another set for test set validation. This approach has more power than the Model_builder_V2, as the full sample is being used to estimate associations. For more information check out: https://www.elderresearch.com/blog/nested-cross-validation-when-cross-validation-isnt-enough/

***
