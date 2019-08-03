# Run_TWAS.R and Run_TWAS.sh

This section contains information on an R script called 'Run_TWAS.R', and a shell script to run the R in an array called 'Run_TWAS.sh'. The R script simply performs a TWAS using FUSION across all chromosomes in parallel and creates a log file.  The shell script just allows you to run the R script for multiple GWAS at a time in an array.

## Pre-requisites
The following software is required for the prediction pipeline:

* GWAS sumstats in LDSC munged format
* FUSION software and reference data (http://gusevlab.org/projects/fusion/)
  * Note this data must be stored in a folder with the following structure:
    * /fusion_twas-master/FUSION.assoc_test.R
    * /SNP-weights/
    * /LDREF/1000G.EUR.
* R packages:
```R
install.packages(c('data.table','foreach','doMC'))
```

## Parameters for Run_TWAS.R
| Flag     | Description                                                  | Default |
| :------- | ------------------------------------------------------------ | :-----: |
| --gwas | GWAS summary statistics in LDSC format [required] | NA |
| --pos | File listing SNP-weights in .pos format [required] | NA |
| --out | Name of output files [required] | NA |
| --fusion_dir | Directory containing fusion software and reference data [required] |   NA    |
| --ncores | Number of cores for parallel analysis [optional] | 5 |

## Parameters for Run_TWAS.sh

| Flag         | Description                                                  | Default |
| :----------- | ------------------------------------------------------------ | :-----: |
| --gwas_list  | A file list GWAS sumstats. First column is GWAS labels, and second column is the path of the summary stats. Columns should be space delimited [required] |   NA    |
| --pos        | File listing SNP-weights in .pos format [required]           |   NA    |
| --outdir     | Name of the directory for output files. A folder will be created for each GWAS [required] |   NA    |
| --fusion_dir | Directory containing fusion software and reference data. [required] |   NA    |
| --ncores     | Number of cores for parallel analysis [required]             |   NA    |

## Output files

A folder will be created for each GWAS/TWAS containing:

* A log file for the Run_TWAS.R script
* A folder of per chromosome log files from FUSION
* Transcriptome-wide results of the TWAS ('_res_GW.txt')
* GWAS summary stats used in the TWAS ('_noNA.sumstats.gz')

## Examples
```sh
# If running for a single GWAS across 8 cores
Rscript Run_TWAS.R \
  --gwas DEPR06.sumstats.gz \
  --pos All_weights.pos \
  --out DEPR06_TWAS \
  --fusion_dir ${fusion_dir} \
  --ncores 8

# If running for 10 GWAS in an array, with 8 cores per GWAS, limiting to 3 at a time.
qsub -l h_vmem=6G -pe smp 8 -tc 3 -t 1-10 Run_TWAS.sh \
  --gwas_list ${TWAS_rep}/todo.txt \
  --pos ${TWAS_rep}/All_tissues.pos \
  --outdir ${TWAS_rep} \
  --fusion_dir ${FUSION_dir} \
  --ncores 8

```
