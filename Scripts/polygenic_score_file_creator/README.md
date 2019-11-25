# polygenic_score_file_creator.R

This is script is used to create reference files for polygenic scoring, including .SCORE files for plink after LD-based clumping, and files containing the mean and standard deviation (SD) for scores for each p-value threshold, across groups within the reference (e.g. ancestries). The script uses PLINK1.9 to perform LD-clumping and scoring in the reference sample. This script was developed for preparing reference data in a genotype-based prediction pipeline. More information [here](https://opain.github.io/GenoPred/Pipeline_prep.html)

## Pre-requisites
The following software is required for the prediction pipeline:

* PLINK v1.9 (https://www.cog-genomics.org/plink2/)
* Per chromosome files for the desired reference genotype data (e.g. 1000 Genomes)
* R packages:
```R
install.packages(c('data.table'))
```

## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --ref_plink_chr | Path to per chromosome reference PLINK files [required] | NA |
| --ref_keep | Keep file to subset individuals in reference for clumping [optional] | NA |
| --sumstats | GWAS summary statistics in LDSC format [required] | 10 |
| --pTs | List of p-value thresholds for scoring [optional] | '1e-8,1e-6,1e-4,1e-2,0.1,0.2,0.3,0.4,0.5,1' |
| --plink | Path PLINK software binary [required] | NA |
| --output | Path for output files [optional] | './Output' |
| --ref_pop_scale | List of keep files for grouping individuals [optional] | NA |
| --memory | Memory limit in Mb [optional] | 5000 |
| --prune_hla | Retain only top associated variant in HLA region [optional] | TRUE |
| --dense | Specify as T for dense thresholding. pTs then interpretted as seq() command wih default 5e-8,1,5e-4 [optional] | TRUE |

## Output files

The script will create per chromosome .score and .range_values files in PLINK format. These files tell PLINK what the effect of each SNP is and the p-value of each SNP to enable selection of SNPs based on p-value thresholds. The script will also create a .NSNP_per_pT files, indicating the number of SNPs surpassing p-value thresholds before and after clumping. A .log file will also be created.

If --ref_pop_scale is specified, the script will also creates files stating the mean and standard deviation of the polygenic scores for each group. 

## Examples
```sh
Rscript polygenic_score_file_creator.R \
  --ref_plink_chr ${Geno_1KG_dir}/1KGPhase3.w_hm3.chr \
  --ref_keep ${Geno_1KG_dir}/keep_files/EUR_samples.keep \
  --sumstats ${gwas_rep}/${gwas}.sumstats.gz \
  --plink ${plink1_9} \
  --memory 3000 \
  --output 1KGPhase3.w_hm3.${gwas} \
  --ref_pop_scale ${Geno_1KG_dir}/super_pop_keep.list
```
