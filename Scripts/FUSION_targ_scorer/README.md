# FUSION_targ_scorer.R

This script is for predicting functional genomic features in a target sample. Instructions for preparing the required reference files for this analysis can be found [here](https://opain.github.io/GenoPred/Pipeline_prep.html#5_functionally-informed_polygenic_scoring). This script is a part of a pipeline for estimating reference standardised genotype-based scores in a target sample, which is further described [here](https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html).

## Pre-requisites
The following software is required for the prediction pipeline:

* PLINK v1.9 (https://www.cog-genomics.org/plink/)
* Per chromosome PLINK files for the target sample that have been previously harmonised with the reference data (more information [here](https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html#1_harmonisation_with_the_reference))
* Reference data for predicting functional genomic features (see [here](https://opain.github.io/GenoPred/Pipeline_prep.html#5_functionally-informed_polygenic_scoring)
* pigz: software for parallel gz compression (https://zlib.net/pigz/)

* R packages:
```R
install.packages(c('data.table','optparse','foreach','doMC'))
```

## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --targ_plink_chr | Path to per chromosome target PLINK binaries [required] | NA |
| --targ_keep | Path to keep file for target [optional] | NA |
| --weights | Path for .pos file describing features [required] | NA |
| --memory | RAM available in MB [required] | 5000 |
| --plink | Path to PLINK software [required] | NA |
| --n_cores | Specify the number of cores available [required] | 1 |
| --score_files | Path to SCORE files corresponding to weights [required] | NA |
| --ref_scale | File containing the scale of features in the reference [required] | NA |
| --ref_freq_chr | Path to per chromosome reference PLINK .frq files [required] | NA |
| --pigz | Path to pigz binary [required] | NA |
| --output | Name of output files [required] | NA |

## Output files
This script will create: 

* A .predictions.gz file containing reference standardised feature predictions
* A .log file

## Examples
```sh
qsub -pe smp 5 Rscript FUSION_targ_scorer.R \
  --targ_plink_chr ${UKBB_output}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr \
  --targ_keep ${UKBB_output}/Projected_PCs/AllAncestry/UKBB.w_hm3.AllAncestry.EUR.keep \
  --weights ${FUSION_dir}/SNP-weights/Blood/Blood.pos \
  --plink ${plink1_9} \
  --n_cores 5 \
  --score_files /users/k1806347/brc_scratch/Data/FUSION/SCORE_FILES \
  --ref_scale ${Geno_1KG_dir}/Predicted_expression/FUSION/Blood/1KGPhase3.w_hm3.FUSION.Blood.EUR.scale \
  --pigz pigz \
  --output ${UKBB_output}/Predicted_expression/FUSION/EUR/Blood/UKBB.w_hm3.QCd.AllSNP.FUSION.Blood \
  --ref_freq_chr ${Geno_1KG_dir}/freq_files/EUR/1KGPhase3.w_hm3.EUR.chr
```
