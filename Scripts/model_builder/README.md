# model_builder

This directory contains the current implementation of GenoPred's predictor evaluation and comparison tooling. It supersedes the earlier `Scripts/Model_builder/Model_builder_V2.R` scripts, using a nested cross-validation design, the `GenoUtils` package, and the shared helper functions in `../functions`.

There are two scripts:

* `model_builder.R` — the main script. Evaluates and compares groups of predictors, building an elastic net model for each group and testing predictive utility with nested cross-validation.
* `model_builder_top1.R` — a variant that focuses on selecting the single best predictor within each group (and optionally a multi-predictor model), using a single layer of cross-validation.

## Overview

`model_builder.R` is used to evaluate and compare predictors (e.g. polygenic scores). Predictors are organised into **groups**, and a model is fitted for each group:

* Groups containing more than one predictor are modelled with an **elastic net** (`glmnet` via `caret`).
* Groups containing a single predictor are modelled with a **GLM**.

Model performance is estimated using **nested cross-validation**: the outer folds provide an unbiased estimate of predictive utility, while the inner folds are used to tune the elastic net hyperparameters. This avoids the optimism that arises when hyperparameters are tuned and evaluated on the same data. The predictive utility of different models is then compared pairwise using the `paired.r` function from the `psych` package, which accounts for the correlation between models' predictions.

The outcome type is detected automatically: a binary outcome is modelled with logistic regression (`binomial`), and an outcome with more than two unique values is modelled with linear regression (`gaussian`).

## Pre-requisites

The scripts source shared helpers from `../functions` using relative paths, so they must be run from within this directory (`Scripts/model_builder`).

Required R packages:

```R
install.packages(c('optparse','data.table','glmnet','doMC','caret','pROC','verification','psych'))
```

The `GenoUtils` package is also required and is installed as part of the GenoPred environment.

## Usage

```bash
Rscript model_builder.R \
  --outcome pheno.txt \
  --predictors predictor_list.txt \
  --out results/my_analysis \
  --n_core 4
```

## Parameters (`model_builder.R`)

| Flag | Description | Default |
| :--- | :--- | :---: |
| --outcome | File containing outcome data [required] | NULL |
| --predictors | File listing files containing predictors, with a group column for model comparison [required] | NULL |
| --n_outer_fold | Number of folds for outer cross-validation (unbiased evaluation) [optional] | 10 |
| --n_inner_fold | Number of folds for inner cross-validation (elastic net tuning) [optional] | 10 |
| --n_core | Number of cores for parallel computing [optional] | 1 |
| --keep | File containing list of individuals to include in analysis [optional] | NULL |
| --outcome_pop_prev | Prevalence of outcome in the general population, used for liability-scale R² [optional] | NULL |
| --out | Prefix for output files [required] | NULL |
| --assoc | Perform association analysis between each predictor and outcome [optional] | TRUE |
| --compare_predictors | Assign each predictor to its own group (in addition to the supplied groups) [optional] | FALSE |
| --pred_miss | Proportion of missing values allowed in a predictor before it is dropped [optional] | 0.1 |
| --top1 | Also evaluate a model using the single best predictor within each group [optional] | FALSE |
| --all_model | Also evaluate a model containing all predictors combined [optional] | TRUE |
| --export_models | Export the coefficients of the final models refitted on all data [optional] | TRUE |
| --seed | Random seed [optional] | 1 |

## Input files

* **`--predictors`**: A table with two columns, `predictor` and `group`.
  * `predictor` — path to a file containing predictor data. Each predictor file should have columns `FID`, `IID`, followed by one or more predictor variables.
  * `group` — the group each predictor file is assigned to. Predictors sharing a group are combined into a single model.
* **`--outcome`**: A file with columns `FID`, `IID`, and a single outcome variable. A binary outcome triggers logistic regression; otherwise linear regression is used.
* **`--keep`** (optional): A file listing `FID` and `IID` of the individuals to retain in the analysis.

Before modelling, predictors are merged across files on `IID`, and predictors are removed if they exceed the `--pred_miss` missingness threshold, have zero variance, or are identical to another predictor within the same group.

## Group construction

Groups are assembled from the `group` column of the predictors file, plus optionally:

* one group per predictor, if `--compare_predictors TRUE`;
* an `all` group containing every predictor, if `--all_model TRUE` and more than one group is present;
* a `top1` model per group, if `--top1 TRUE`, which selects the single best predictor within a group using the training data and evaluates it in the held-out data.

## Output files

All outputs are prefixed with the value supplied to `--out`.

* **`.group_list.txt`** — the groups analysed and the number of predictors in each.
* **`.assoc.txt`** — univariate association of each predictor with the outcome (written when `--assoc TRUE`). Columns:
  * `Group` — the group the predictor belongs to
  * `Predictor` — predictor name
  * `BETA` — effect size, in standard-deviation units
  * `SE` — standard error of `BETA`
  * `P` — association p-value
  * `Obs_R2` — variance explained on the observed scale
  * `N` — sample size
  * For binary outcomes, additionally: `N_case`, `N_control`, and `Liab_R2` (variance explained on the liability scale, derived from `--outcome_pop_prev`).
* **`.pred_eval.txt`** — predictive utility of each model, estimated from the pooled out-of-fold predictions (`Group` plus evaluation metrics).
* **`.pred_comp.txt`** — pairwise comparison of models (written when more than one group is present). Columns:
  * `Model_1`, `Model_2` — the two models compared
  * `Model_1_R`, `Model_2_R` — correlation between each model's predictions and the outcome
  * `R_diff` — difference in correlation
  * `R_diff_pval` — p-value for the difference, from `paired.r`
* **`.log`** — run log.
* **`final_models/`** — coefficients of the final models refitted on all data (written when `--export_models TRUE`), saved in the same directory as `--out`.

## `model_builder_top1.R`

This variant is used when the aim is to identify the single best predictor within each group, optionally alongside a combined multi-predictor model. It uses a single layer of cross-validation rather than the nested design.

The predictors file for this script uses columns `predictor`, `multi`, and `top1`, where `multi` and `top1` flag whether each predictor file participates in multi-predictor modelling and/or best-predictor selection.

### Parameters (`model_builder_top1.R`)

| Flag | Description | Default |
| :--- | :--- | :---: |
| --outcome | File containing outcome data [required] | NULL |
| --predictors | File listing predictor files, with `multi` and `top1` columns [required] | NULL |
| --n_fold | Number of folds for cross-validation [optional] | 10 |
| --n_core | Number of cores for parallel computing [optional] | 1 |
| --keep | File containing list of individuals to include in analysis [optional] | NULL |
| --outcome_pop_prev | Prevalence of outcome in the general population [optional] | NULL |
| --out | Prefix for output files [required] | NULL |
| --pred_miss | Proportion of missing values allowed in a predictor [optional] | 0.1 |
| --export_models | Export final model coefficients [optional] | TRUE |
| --seed | Random seed [optional] | 1 |

Exported final models are written to `<out>_final_models`.

---

*Scripts written by Oliver Pain. This README documents the current `Scripts/model_builder` implementation; the earlier single-split implementation is retained under `Scripts/Model_builder`.*
