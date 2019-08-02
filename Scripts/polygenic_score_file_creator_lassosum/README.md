# polygenic_score_file_creator_lassosum.R

This is script is used to create reference files for polygenic scoring using lassosum. For more inf. This script was developed for preparing reference data in a genotype-based prediction pipeline. More information on lassosum can be found [here](https://github.com/tshmak/lassosum).

## Pre-requisites
The following software is required for the prediction pipeline:

* PLINK v1.9 (https://www.cog-genomics.org/plink2/)
* Per chromosome files for the desired reference genotype data (e.g. 1000 Genomes)
* lassosum (https://github.com/tshmak/lassosum)
* R packages:
```R
install.packages(c('data.table'))
```

## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --ref_plink_gw | Path to genome-wide reference PLINK files [required] | NA |
| --ref_keep | Keep file to subset individuals in reference [optional] | NA |
| --sumstats | GWAS summary statistics in LDSC format [required] | 10 |
| --output | Path for output files [optional] | './Output' |
| --ref_pop_scale | List of keep files for grouping individuals [optional] | NA |

## Output files

The script will create .unvalidated.model.RDS file containing a lassosum object for scoring across a range of shrinkage parameters. The script will also create a .pseudovalidated.model.RDS file containing a lassosum object for scoring based on a single set of shrinkage parameters, predicted to be the best fit according to pseudovalidation.

If --ref_pop_scale is specified, the script will also creates files stating the mean and standard deviation of the polygenic scores for each group. 

## Examples
```sh
Rscript polygenic_score_file_creator_lassosum.R \
	--ref_plink_gw ${Geno_1KG_dir}/1KGPhase3.w_hm3.GW \
	--ref_keep ${Geno_1KG_dir}/keep_files/EUR_samples.keep \
	--sumstats ${gwas_rep}/${gwas}.sumstats.gz \
	--output 1KGPhase3.w_hm3.${gwas} \
	--ref_pop_scale ${Geno_1KG_dir}/super_pop_keep.list
```
