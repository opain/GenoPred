# Ancestry_identifier.R

This is script is used to infer the global ancestry of a target sample based on a reference-projected principal component (PC) elastic net model.

## Pre-requisites
The following software is required for the prediction pipeline:

* PLINK v1.9 (https://www.cog-genomics.org/plink2/).
* PLINK v2 (https://www.cog-genomics.org/plink/2.0/).
* Per chromosome files for the desired reference genotype data in binary PLINK format (e.g. 1000 Genomes).
* Per chromsome or genome-wide target genotype data in binary PLINK format.
* Per chromsome or genome-wide target genotype data in binary PLINK format.
* A file listing the path to a series of keep files for the reference genotype data used for scaling PCs.
* A file indicating which population reference individuals correspond to.

See example code below and corresponding files to see examples.

* R packages:
```R
install.packages(c('data.table','caret','pROC','verification','ggplot2','cowplot'))
```

## Parameters

| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --target_plink_chr | Path to per chromosome target PLINK files [required] | NA |
| --target_plink | Path to genome-wide target PLINK files [required] | NA |
| --ref_plink_chr | Path to per chromosome reference PLINK files [required] | NA |
| --target_fam | Target sample fam file. [optional] | NA |
| --maf | Minor allele frequency threshold [optional] | NA |
| --geno | Variant missingness threshold [optional] | 0.02 |
| --hwe | Hardy Weinberg p-value threshold. [optional] | NA |
| --n_pcs | Number of PCs (min=4) [optional] | 10 |
| --plink | Path PLINK software binary [required] | plink |
| --plink2 | Path PLINK software binary [required] | plink |
| --output | Path for output files [required] | ./PC_projector_output/Output |
| --ref_pop_scale | List of keep files for ancestry specific scaling [optional] | NA |
| --pop_data | Population data for the reference samples [optional] | NA |
| --model_method | Method used for generate prediction model. Only glmnet tested. [optional] | glmnet |
| --SD_rule | Logical indicating whether the 3SD rule should be used to define ancestry, or the model-based approach [optional] | FALSE |
| --prob_thresh | Indicates whether probability threshold should be used when defining ancestry [optional] | NA |
| --memory | Memory limit [optional] | 5000 |


## Output files

The script will calculate PCs for reference individuals (output.eigenvec), the SNP weights for the PCs (output.eigenvec.var), the mean and SD of each PC in the full reference (output.scale), and the mean and SD of each PC in each reference group specified in ref_pop_scale (output.group.scale). The script will also output the elastic net model predicting groups in pop_data as an .rds, along with some accuracy metrics (output.pop_data_prediction _details.txt). 

The script will also files indicating the ancestry of the target sample. The script will output keep files indicating the ancestry group of each individual, determined using the 3SD rule (output.ref_pop_scale.keep), and the PC scores scaled and centred using the reference mean and SD (output.pop_data.eigenvec). The script will also output keep files indicating the ancestry of the target samples based on the elastic net model (output.model_pred.pop_data). The script will output .png files showing the PCs of target samples relative to the reference samples.

## Examples
```sh
Rscript Ancestry_identifier.R \
  --target_plink_chr target/target.chr \
  --ref_plink_chr ref/ref.chr \
  --n_pcs 6 \
  --plink plink1.9 \
  --plink2 plink2 \
  --output test \
  --ref_pop_scale ref_super_pop_keep_list.txt \
  --pop_data ref_pop_dat.txt \
  --prob_thresh 0.5
```
