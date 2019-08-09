# Scaled_polygenic_scorer_lassosum.R

This script is for calculating reference standardised polygenic scores in a target sample. Polygenic scores are calculated using the [lassosum](https://github.com/tshmak/lassosum). Instructions for preparing the required reference files for this analysis can be found [here](https://opain.github.io/GenoPred/Pipeline_prep.html#43_prepare_score_files_and_scaling_files_for_polygenic_scoring_using_lassosum). This script is a part of a pipeline for estimating reference standardised genotype-based scores in a target sample, which is further described [here](https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html).

## Pre-requisites
The following software is required for the prediction pipeline:

* lassosum (https://github.com/tshmak/lassosum)
* Per chromosome PLINK files for the target sample that have been previously harmonised with the reference data (more information [here](https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html#1_harmonisation_with_the_reference))
* Reference data for polygenic scoring (lassosum method) (see [here](https://opain.github.io/GenoPred/Pipeline_prep.html#43_prepare_score_files_and_scaling_files_for_polygenic_scoring_using_lassosum))

* R packages:
```R
install.packages(c('data.table','optparse','lassosum','foreach','doMC'))
```

## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --target_plink_chr | Path to per chromosome target PLINK files [required] | NA |
| --target_keep | Path to keep file for target [optional] | NA |
| --ref_model | Path to reference scoring files [required] | NA |
| --output | Path for output files [required] | ./Output |
| --ref_scale | Path reference scale file [required] | NA |
| --n_cores | Path reference scale file [required] | 1 |
| --pheno_name | Name of phenotype to be added to column names. Default is SCORE. [optional] | ./Output |

## Output files
This script will create: 

* A .lassosum_profiles file containing reference standardised polygenic scores
* A .log file 

## Examples
```sh
qsub -pe smp 6 Rscript Scaled_polygenic_scorer_lassosum.R \
  --target_plink_chr ${UKBB_output}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr \
  --target_keep ${UKBB_output}/Projected_PCs/AllAncestry/UKBB.w_hm3.AllAncestry.EUR.keep \
  --ref_model ${Geno_1KG_dir}/Score_files_for_poylygenic_lassosum/DEPR06/1KGPhase3.w_hm3.DEPR06.unvalidated.model.RDS \
  --ref_scale ${Geno_1KG_dir}/Score_files_for_poylygenic_lassosum/DEPR06/1KGPhase3.w_hm3.DEPR06.EUR.scale \
  --pheno_name DEPR06 \
  --n_cores 6 \
  --output ${UKBB_output}/PolygenicScores_lassosum/EUR/DEPR06/UKBB.w_hm3.DEPR06.EUR
```
