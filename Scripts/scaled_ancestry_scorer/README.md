# scaled_ancestry_scorer.R

This script can be used to project principal components of population structure (mainly ancestry) that have been previously identified in another sample (such as the 1000 Genomes reference). The script scales the projected PCs according to the reference, can assign individuals to groups in the reference (ancestral groups), and create plots to compare the target and reference individuals.  Instructions for preparing the required reference files for this analysis can be found [here](https://opain.github.io/GenoPred/Pipeline_prep.html#3_ancestry_scoring). This script is a part of a pipeline for estimating reference standardised genotype-based scores in a target sample, which is further described [here](https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html).

## Pre-requisites
The following software is required for the prediction pipeline:

* PLINK v2 (https://www.cog-genomics.org/plink/2.0/)
* Per chromosome PLINK files for the target sample that have been previously harmonsied with the reference data (more information [here](https://opain.github.io/GenoPred/Genotype-based_scoring_in_target_samples.html#1_harmonisation_with_the_reference))
* Reference data for principal component projection (see [here](https://opain.github.io/GenoPred/Pipeline_prep.html))

* R packages:
```R
install.packages(c('data.table','optparse','ggplot2','cowplot'))
```

## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --target_plink_chr | Path to per chromosome target PLINK files [required] | NA |
| --target_keep | Path to keep file for target sample individuals [optional] | NA |
| --ref_freq_chr | Path to per chromosome reference PLINK .frq files [required] | NA |
| --ref_eigenvec | Path to reference PLINK .eigenvec file [required] | NA |
| --ref_score | Path to per chromosome reference PLINK .eigenvec.var file [required] | NA |
| --plink2 | Path PLINK software binary [required] | plink |
| --output | Path for output files [required] | ./PC_projector_output/Output |
| --memory | Memory limit [optional] | 5000 |
| --pop_data | Path to file containing reference population data [optional] | NA |
| --pop_model | rds file containing population prediction model [optional] | NA |
| --pop_model_scale | Path to reference scaling file used when deriving the pop_model [optional] | NA |
| --pop_scale_for_keep | Path to population specific scale files for excluding outliers [optional] | NA |
| --ref_scale | Path to population specific scale file [optional] | NA |

Note. Either --ref_scale or --pop_scale_for_keep must be specified to scale the eigenvectors.

## Output files
This script will always create: 

* .eigenvec file containing principal component scores for individuals in the target sample for each PC in the reference data. 
* .log file

If --pop_data is specified, the script will also create figures comparing the scores of the target sample, to those in the reference sample. 
If --pop_model is specified, the script will create a .model_pred file containing the probability of each individual belong to each group in the model.
If --pop_scale_for_keep is specified, the script will create a .keep file listing individuals within ±3SD of the mean of each group specified, and create an .eigenvec file for each group.

## Examples
```sh
####
# Calculate all ancestry PCs of population structure and predict ancestry
####

qsub Rscript scaled_ancestry_scorer.R \
	--target_plink_chr ${UKBB_output}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr \
	--ref_freq_chr ${Geno_1KG_dir}/freq_files/AllAncestry/1KGPhase3.w_hm3.AllAncestry.chr \
	--ref_score ${Geno_1KG_dir}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.eigenvec.var \
	--ref_eigenvec ${Geno_1KG_dir}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.eigenvec \
	--plink2 ${plink2} \
	--output ${UKBB_output}/Projected_PCs/AllAncestry/UKBB.w_hm3.AllAncestry \
	--pop_data ${Geno_1KG_dir}/integrated_call_samples_v3.20130502.ALL.panel_small \
	--pop_model ${Geno_1KG_dir}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.pop_enet_model.rds \
	--pop_model_scale ${Geno_1KG_dir}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.scale \
	--pop_scale_for_keep ${Geno_1KG_dir}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry.scalefiles.txt

####
# Calculate European ancestry PCs of population structure
####

# Use the EUR.keep file from the AllAncestry run above to select people of EUR ancestry.
qsub Rscript scaled_ancestry_scorer.R \
	--target_plink_chr ${UKBB_output}/Genotype/Harmonised/UKBB.w_hm3.QCd.AllSNP.chr \
	--target_keep ${UKBB_output}/Projected_PCs/AllAncestry/UKBB.w_hm3.AllAncestry.EUR.keep \
	--ref_freq_chr ${Geno_1KG_dir}/freq_files/EUR/1KGPhase3.w_hm3.EUR.chr \
	--ref_score ${Geno_1KG_dir}/Score_files_for_ancestry/EUR/1KGPhase3.w_hm3.EUR.eigenvec.var \
	--ref_eigenvec ${Geno_1KG_dir}/Score_files_for_ancestry/EUR/1KGPhase3.w_hm3.EUR.eigenvec \
	--ref_scale ${Geno_1KG_dir}/Score_files_for_ancestry/EUR/1KGPhase3.w_hm3.EUR.scale \
	--plink2 ${plink2} \
	--output ${UKBB_output}/Projected_PCs/EUR/UKBB.w_hm3 \
	--pop_data ${Geno_1KG_dir}/integrated_call_samples_v3.20130502.ALL.panel_small
```
