# Harmonisation_of_UKBB.R

This script is used to harmonise the UK-Biobank bank genetic data with the reference and format it to be ready for genotype-based scoring. This script was developed for preparing target data in a genotype-based prediction pipeline. More information [here](https://opain.github.io/GenoPred/Pipeline_prep.html) 

## Pre-requisites
The following software is required for the prediction pipeline:

* PLINK v1.9 (https://www.cog-genomics.org/plink2/)

* QCTOOL V2 (https://www.well.ox.ac.uk/~gav/qctool_v2)

* Per chromosome imputed UK-biobank genetic data in .bgen format (as received from UK-biobank)

* R packages:
```R
install.packages(c('data.table','optparse'))
```

## Parameters
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --chr | Chromosome number [required] | NA |
| --input_dir | Directory containing imputed UKBB data [required] | NA |
| --target_fam | Path to fam file for target sample [required] | NA |
| --reference_dir | Directory containing reference data [required] | NA |
| --plink | Path to PLINK V1.9 software binary [required] | NA |
| --qctool2 | Path to QCTOOL V2 binary [required] | NA |
| --output_dir | Output folder name [required] | NA |
| --debug | Set to T to create object for debugging [optional] | F |

## Output files

The script will create a set of PLINK .bed/.bim/.fam files in the output_dir for each chromosome, containing only SNPs within the HapMap3 SNP list, flipped to match the strand in the reference, and with missing SNPs inserted as NA.

## Examples
```sh
for chr in $(seq 1 22); do
qsub -l h_vmem=10G Rscript Harmonisation_of_UKBB.R \
  --chr $chr \
  --input_dir UKBB/Imputed \
  --target_fam UKBB/Imputed/UKBB.sample \
  --output_dir UKBB/Harmonised \
  --reference_dir 1KG/Phase3 \
  --qctool2 qctool2 \
  --plink plink
done
```
