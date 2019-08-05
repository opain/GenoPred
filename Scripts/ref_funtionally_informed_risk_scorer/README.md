# polygenic_score_file_creator.R

This is script is used to create reference files for functionally informed polygenic scoring based on predicted features (e.g. gene expression) and inferred feature association summary statistics (e.g. TWAS). This script was developed for preparing reference data in a genotype-based prediction pipeline. More information [here](https://opain.github.io/GenoPred/Pipeline_prep.html)

## Pre-requisites
The following software is required for the prediction pipeline:

* Feature predictions in a reference sample (See [here](https://github.com/opain/GenoPred/tree/master/Scripts/FUSION_ref_scorer))
* Feature associations (e.g. TWAS) results from FUSION (http://gusevlab.org/projects/fusion/)
* R packages:
```R
install.packages(c('optparse','data.table'))
```
## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --twas_results | File containing TWAS results [required] | NA |
| --ref_feature_pred | File containing feature predictions for the reference sample [required] | NA |
| --output | Name of output files [required] | NA |
| --clump_thresh | r2 threshold for clumping [optional] | 0.9 |
| --cor_window | Window for clumping [optional] | 5e6 |
| --pTs | Window for clumping [optional] | '1,5e-1,1e-1,5e-2,1e-2,1e-3,1e-4,1e-5,1e-6' |
| --clump_mhc | Retain only the most significant gene within the MHC region [optional] | TRUE |
| --ref_keep | Keep file for reference individuals [optional] | NA |
| --panel | Panel from TWAS [optional] | NA |
| --r2_weighted | Set to T if gene expression should be weighted by R2 of predicted expression [optional] | F |
| --ref_scale | Path to file for scaling feature predictions [required] | NA |

## Output files

The script will create:
* .score file containing effect size and p value for each feature after LD-based clumping
* .scale file containing the mean and SD of scores within the reference sample
* .NFeat file showing the number of features included at different p-value thresholds
* .log file

## Examples
```sh
Rscript ref_funtionally_informed_risk_scorer.R \
  --twas_results DEPR06_res_GW.txt \
  --ref_feature_pred 1KGPhase3.w_hm3.FUSION.Blood.predictions.gz \
  --output DEPR06/1KGPhase3.w_hm3.EUR.FUSION.DEPR06.Blood \
  --ref_keep EUR_samples.keep \
  --ref_scale Predicted_expression/FUSION/Blood/1KGPhase3.w_hm3.FUSION.Blood.EUR.scale \
  --panel Blood
```
