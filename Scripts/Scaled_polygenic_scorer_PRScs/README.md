# Scaled_polygenic_scorer_PRScs.R

This script is for calculating reference standardised polygenic scores in a target sample. Polygenic scores are calculated using [PRScs](https://github.com/getian107/PRScs). Instructions for preparing the required reference files for this analysis can be found [here](https://opain.github.io/GenoPred/Pipeline_prep.html#44_prepare_score_and_scale_files_for_polygenic_scoring_using_prscs). This script is a part of a pipeline for estimating reference standardised genotype-based scores in a target sample, which is further described [here](https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html).

## Pre-requisites
The following software is required for the prediction pipeline:

* PLINK v1.9 (https://www.cog-genomics.org/plink/)
* Per chromosome PLINK files for the target sample that have been previously harmonised with the reference data (more information [here](https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html#1_harmonisation_with_the_reference))
* Reference data for polygenic scoring (pT + clump method) (see [here](https://opain.github.io/GenoPred/Pipeline_prep.html#44_prepare_score_and_scale_files_for_polygenic_scoring_using_prscs))

* R packages:
```R
install.packages(c('data.table','optparse'))
```

## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --target_plink_chr | Path to per chromosome target PLINK files [required] | NA |
| --target_keep | Path to keep file for target [optional] | NA |
| --ref_score | Path to reference scoring files [required] | NA |
| --ref_freq_chr | Path to per chromosome reference PLINK .frq files [required] | NA |
| --plink | Path PLINK v1.9 software binary [required] | plink |
| --output | Path for output files [required] | ./Output |
| --ref_scale | Path reference scale file [required] | NA |
| --pheno_name | Name of phenotype to be added to column names. Default is SCORE. [optional] | ./Output |
| --memory | Memory limit [optional] | 5000 |

## Output files
This script will create: 

* A .PRScs_profiles file containing reference standardised polygenic scores
* A .log file 

## Examples
```sh
qsub Rscript Scaled_polygenic_scorer_PRScs/Scaled_polygenic_scorer_PRScs.R \
  --target_plink_chr ${UKBB_output}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr \
  --target_keep ${UKBB_output}/Projected_PCs/AllAncestry/UKBB.w_hm3.AllAncestry.EUR.keep \
  --ref_score ${Geno_1KG_dir}/Score_files_for_poylygenic_PRScs/DEPR06/1KGPhase3.w_hm3.DEPR06 \
  --ref_scale ${Geno_1KG_dir}/Score_files_for_poylygenic_PRScs/DEPR06/1KGPhase3.w_hm3.DEPR06.EUR.scale \
  --ref_freq_chr ${Geno_1KG_dir}/freq_files/EUR/1KGPhase3.w_hm3.EUR.chr \
  --plink plink \
  --pheno_name DEPR06 \
  --output ${UKBB_output}/PolygenicScores_PRScs/EUR/DEPR06/UKBB.w_hm3.DEPR06.EUR
```
