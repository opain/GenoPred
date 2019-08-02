# polygenic_score_file_creator.R

This is script is used to create reference files for polygenic scoring based on GWAS sumstats processed by PRScs, a method for performing Bayesian continuous shrinkage to GWAS sumstats to account for LD and winners curse. More information on PRScs can be found [here](https://github.com/getian107/PRScs) . This script was developed for preparing reference data in a genotype-based prediction pipeline. More information [here](https://opain.github.io/GenoPred/Pipeline_prep.html)

## Pre-requisites
The following software is required for the prediction pipeline:

* PLINK v1.9 (https://www.cog-genomics.org/plink2/)
* Per chromosome files for the desired reference genotype data (e.g. 1000 Genomes)
* PRScs software and reference data (https://github.com/getian107/PRScs)
* Python V2.X
* R packages:
```R
install.packages(c('data.table','foreach','doMC'))
```

## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --ref_plink_chr | Path to per chromosome reference PLINK files [required] | NA |
| --ref_keep | Keep file to subset individuals in reference for clumping [optional] | NA |
| --sumstats | GWAS summary statistics in LDSC format [required] | 10 |
| --PRScs_path | List of p-value thresholds for scoring [required] | NA |
| --PRScs_ref_path | Path to PRScs reference [required] | NA |
| --phi_param | Path to PRScs reference [optional] | 'auto' |
| --plink | Path PLINK software binary [required] | NA |
| --output | Path for output files [optional] | './Output' |
| --ref_pop_scale | List of keep files for grouping individuals [optional] | NA |
| --memory | Memory limit in Mb [optional] | 5000 |
| --prune_hla | Retain only top associated variant in HLA region [optional] | TRUE |
| --n_cores | Number of cores for parallel computing [optional] | 1 |
| --python_path | Path to python 2.X [required] | NA |

## Output files

The script will create per chromosome and per phi_param .score files in PLINK format. A .log file will also be created.

If --ref_pop_scale is specified, the script will also creates files stating the mean and standard deviation of the polygenic scores for each group. 

## Examples
```sh
Rscript.sh polygenic_score_file_creator_PRScs.R \
--ref_plink_chr ${Geno_1KG_dir}/1KGPhase3.w_hm3.chr \
--ref_keep ${Geno_1KG_dir}/keep_files/EUR_samples.keep \
--sumstats ${gwas_rep}/${gwas}.sumstats.gz \
--plink ${plink1_9} \
--memory 5000 \
--output 1KGPhase3.w_hm3.${gwas} \
--ref_pop_scale ${Geno_1KG_dir}/super_pop_keep.list \
--PRScs_path ${PRScs_dir}/PRScs.py \
--PRScs_ref_path ${PRScs_dir}/ldblk_1kg_eur \
--n_cores 6 \
--phi_param 1e-6,1e-4,1e-2,1,auto
```
