# polygenic_score_file_creator.R

This is script is used to create reference files for polygenic scoring, including .SCORE files for plink after LD-based clumping, and files containing the mean and standard deviation (SD) for scores for each p-value threshold, across groups within the reference (e.g. ancestries). The script uses PLINK1.9 to perform LD-clumping and scoring in the reference sample.



This script is used to generate principal components (PCs) of ancestry, and the files required for projecting the derived PCs onto target samples. The script uses PLINK1.9 and PLINK2. This script was developed for preparing reference data in a genotype-based prediction pipeline. More information [here](https://opain.github.io/GenoPred/Pipeline_prep.html)

## Pre-requisites
The following software is required for the prediction pipeline:

* PLINK v1.9 (https://www.cog-genomics.org/plink2/)

* PLINK v2 (https://www.cog-genomics.org/plink/2.0/)

* Per chromosome files for the desired reference genotype data (e.g. 1000 Genomes)

* R packages:
```R
install.packages(c('data.table','caret','pROC','verification'))
```

## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --ref_plink_chr | Path to per chromosome reference PLINK files [required] | NA |
| --ref_keep | Keep file to subset individuals in reference for clumping [optional] | NA |
| --n_pcs | Number of PCs [optional] | 10 |
| --plink | Path PLINK software binary [required] | NA |
| --plink2 | Path PLINK2 software binary [required] | NA |
| --output | Path for output files [optional] | './PC_projector_output/Output' |
| --ref_pop_scale | List of keep files for grouping individuals [optional] | NA |
| --memory | Memory limit in Mb [optional] | 5000 |

## Output files

The script will always create .eigenvec and .eigenvec.var files. These contain the PCs score for individuals in the reference dataset, and the SNP-weights for the PCs respectively. 

If --ref_pop_scale is specified, the script will also creates files stating the mean and standard deviation of the PCs for each group. Furthermore, it will derive an elastic net model predicting each group, and report the accuracy of the derived models.

## Examples
```sh
# Create score files for European specific PCs, and scaling files for Europeans
qsub /users/k1806347/brc_scratch/Software/Rscript.sh /users/k1806347/brc_scratch/Software/MyGit/GenoPred/Scripts/ancestry_score_file_creator/ancestry_score_file_creator.R \
	--ref_plink_chr ${Geno_1KG_dir}/1KGPhase3.w_hm3.chr \
	--ref_keep ${Geno_1KG_dir}/keep_files/EUR_samples.keep \
	--plink ${plink1_9} \
	--plink2 ${plink2} \
	--n_pcs 100 \
	--output ${Geno_1KG_dir}/Score_files_for_ancestry/EUR/1KGPhase3.w_hm3.EUR

# Create score files for all ancestry PCs, and scaling files for each super population
qsub /users/k1806347/brc_scratch/Software/Rscript.sh /users/k1806347/brc_scratch/Software/MyGit/GenoPred/Scripts/Pipeline_prep/ancestry_score_file_creator.R \
	--ref_plink_chr ${Geno_1KG_dir}/1KGPhase3.w_hm3.chr \
	--plink ${plink1_9} \
	--plink2 ${plink2} \
	--n_pcs 100 \
	--ref_pop_scale ${Geno_1KG_dir}/super_pop_keep.list \
	--output ${Geno_1KG_dir}/Score_files_for_ancestry/AllAncestry/1KGPhase3.w_hm3.AllAncestry
```
