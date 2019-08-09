# Scaled_polygenic_scorer.R

This script is for calculating reference standardised polygenic scores in a target sample. Polygenic scores are calculated using the traditional p-value thresholding and LD-based clumping approach. Instructions for preparing the required reference files for this analysis can be found [here](https://opain.github.io/GenoPred/Pipeline_prep.html#42_prepare_score_files_and_scaling_files_for_polygenic_scoring_(pt_+_clump)). This script is a part of a pipeline for estimating reference standardised genotype-based scores in a target sample, which is further described [here](https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html).

## Pre-requisites
The following software is required for the prediction pipeline:

* Feature predictions files for the target sample. (more information [here](https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html#41_predicting_functional_genomic_features)
* Reference data for functionally informed polygenic scoring (see [here](https://opain.github.io/GenoPred/Pipeline_prep.html#5_functionally-informed_polygenic_scoring))

* R packages:
```R
install.packages(c('data.table','optparse','foreach','doMC'))
```

## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --targ_feature_pred | Path to predicted expression file [required] | NA |
| --ref_score | Path to reference scoring files [required] | NA |
| --output | Path for output files [required] | ./Output |
| --ref_scale | Path reference scale file [required] | NA |
| --pheno_name | Name of phenotype to be added to column names. Default is SCORE.[optional] | NA |
| --pigz | Path to pigz binary [required] | NA |
| --batch_size | Number of individuals to be processed in each batch [optional] | 5000 |
| --n_cores | Specify the number of cores available [optional] | 1 |

## Output files
This script will create: 

* A .fiprofile file containing reference-standardised functionally-informed polygenic scores
* A .log file 

## Examples
```sh
qsub -pe smp 5 Rscript scaled_functionally_informed_risk_scorer.R \
  --targ_feature_pred ${UKBB_output}/Predicted_expression/FUSION/EUR/Blood/UKBB.w_hm3.QCd.AllSNP.FUSION.Blood.predictions.gz \
  --ref_score ${Geno_1KG_dir}/Score_files_for_functionally_informed_risk_scores/DEPR06/1KGPhase3.w_hm3.EUR.FUSION.DEPR06.Blood.score \
  --ref_scale ${Geno_1KG_dir}/Score_files_for_functionally_informed_risk_scores/DEPR06/1KGPhase3.w_hm3.EUR.FUSION.DEPR06.Blood.scale \
  --pheno_name DEPR06 \
  --n_cores 5 \
  --pigz pigz \
  --output ${UKBB_output}/FunctionallyInformedPolygenicScores/EUR/Blood/UKBB.w_hm3.EUR.Blood.DEPR06
```
