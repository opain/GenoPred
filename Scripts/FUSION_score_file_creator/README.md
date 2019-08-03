# FUSION_score_file_creator.R

This script is to convert FUSION SNP-weights into .score files in PLINK format for predicting features.

## Pre-requisites
The following software is required for the prediction pipeline:

* FUSION SNP-weights (http://gusevlab.org/projects/fusion/)
* R packages:
```R
install.packages(c('data.table','foreach','doMC'))
```

## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --weights | Path for .pos file describing features [required] | NA |
| --weights_dir | Directory containing the weights corresponding to the features in the .pos file [required] | NA |
| --n_cores | Specify the number of cores available [optional]             | 1 |
| --output | Name of output directory [required] | NA |

## Output files

The script will create a folder containing a score file and a snplist file for every feature in the .pos file.

## Examples
```sh
qsub -pe smp 15 -l h_vmem=2G Rscript FUSION_score_file_creator.R \
  --weights All_tissues.pos \
  --weights_dir SNP-weights \
  --output SCORE_FILES \
  --n_cores 15
```
