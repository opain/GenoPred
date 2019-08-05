# FUSION_ref_scorer.R

This script is for predicting functional genomic features (e.g. gene expression) based on genotypic data and FUSION SNP-weights. This script was developed for preparing reference data in a genotype-based prediction pipeline. More information [here](https://opain.github.io/GenoPred/Pipeline_prep.html)

## Pre-requisites
The following software is required for the prediction pipeline:

* PLINK v1.9 (https://www.cog-genomics.org/plink2/)
* Per chromosome files for the desired reference genotype data (e.g. 1000 Genomes)
* FUSION SNP-weights in .score format (See [here](https://github.com/opain/https://github.com/opain/GenoPred/tree/master/Scripts/FUSION_score_file_creator))
* GWAS summary statistics in LDSC munged format
* pigz: software for parallel gz compression (https://zlib.net/pigz/)
* R packages:
```R
install.packages(c('data.table','foreach','doMC'))
```

make_option("--ref_plink_chr", action="store", default=NA, type='character',
		help="Path to per chromosome reference PLINK binaries [required]"),
make_option("--weights", action="store", default=NA, type='character',
		help="Path for .pos file describing features [required]"),
make_option("--n_cores", action="store", default=1, type='numeric',
		help="Specify the number of cores available [optional]"),
make_option("--score_files", action="store", default=NA, type='character',
		help="Path to SCORE files corresponding to weights [required]"),
make_option("--memory", action="store", default=5000, type='numeric',
		help="RAM available in MB [required]"),
make_option("--plink", action="store", default='NA', type='character',
		help="Path to PLINK software [required]"),
make_option("--output", action="store", default=NA, type='character',
		help="Name of output directory [required]"),
make_option("--ref_pop_scale", action="store", default=NA, type='character',
		help="File containing the population code and location of the keep file [required]"),
make_option("--pigz", action="store", default=NA, type='character',
		help="Path to pigz binary [required]")
)

## Parameters
| Flag     | Description                                                  | Default |
| :------------- | ------------------------------------------------------ | :----------: |
| --ref_plink_chr | Path to per chromosome reference PLINK files [required] | NA |
| --plink | Path PLINK software binary [required] | NA |
| --output | Path for output files [required] | NA |
| --ref_pop_scale | List of keep files for grouping individuals [optional] | NA |
| --memory | Memory limit in Mb [optional] | 5000 |
| --weights | Path for .pos file describing features [required] | NA |
| --n_cores | Specify the number of cores available [optional] | 1 |
| --score_files | Path to SCORE files corresponding to weights [required] | NA |
| --pigz | Path to pigz binary [required] | NA |

## Output files

The script will create a .predictions.gz file containing feature predictions, a Prediction_failed.txt file listing features that could not be predicted and the reason, and a .log file.

If --ref_pop_scale is specified, the script will also creates files stating the mean and standard deviation of each feature for each group.

## Examples
```sh
Rscript FUSION_ref_scorer.R \
  --ref_plink_chr 1KGPhase3.w_hm3.chr \
  --weights Brain.pos \
  --output ${Geno_1KG_dir}/Predicted_expression/FUSION/${weights}/1KGPhase3.w_hm3.FUSION.Brain \
  --plink plink1.9 \
  --n_cores 5 \
  --pigz ~/Software/pigz \
  --score_files FUSION/SCORE_FILES \
  --ref_pop_scale super_pop_keep.list
```
