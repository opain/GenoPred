# Genotype-based Prediction Pipeline (GenoPredPipe)

This is a snakemake pipeline for running the GenoPred scripts. The pipeline currently identifies the ancestry of each individual, calculates ancestry-matched reference-projected principal components, and calculates polygenic scores using a range of methods standardised using an ancestry matched reference. Finally, the pipeline provides a report of the results either for each individual or the sample overall.

***

## Getting started

### Step 1

Clone the repository

```bash
git clone git@github.com:https://github.com/opain/GenoPred.git
cd GenoPred/GenoPredPipe
```

### Step 2

Install [Anaconda](https://conda.io/en/latest/miniconda.html).

**Linux:**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
```

Install Python 3.8, Snakemake 5.32, and the basic project dependencies.

```bash
conda activate base
conda install python=3.8
conda install -c conda-forge mamba
mamba install -c bioconda -c conda-forge snakemake-minimal==5.32.2
```

### Step 3

Prepare a snakemake profile for parallel computing. Follow instructions [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles). 

***

## Running pipeline

### Example using test data

#### Step 1

Download test data.

```bash
wget -O test_data.tar.gz https://zenodo.org/record/5205711/files/test_data.tar.gz?download=1
tar -xf test_data.tar.gz
rm test_data.tar.gz
```

#### Step 2

Run full pipeline.

```bash
snakemake --profile slurm --use-conda run_create_reports
```

***

### Using your own data

You must specify a file listing target samples (target_list) and a file listing of GWAS summary statistics (gwas_list) for the pipeline to use. The location of these files should be specified in the config.yaml file. The column names for the gwas_list and target_list files should be as follows:

#### target_list

- name: Name of the target sample
- path: File path to the target genotype data
- type: Either '23andMe' for 23andMe formatted data for an individual, or 'samp_imp' for PLINK1 binary format data (.bed/.bim/.fam) for a sample that has already undergone genotype imputation.

#### gwas_list

- name: Short name for the GWAS
- path: File path to the GWAS summary statistics (uncompressed or gzipped)
- population: The super population that the GWAS was performed in (AFR/AMR/EAS/EUR/SAS)
- sampling: The proportion of the GWAS sample that were cases (if binary, otherwise NA)
- prevalence: The population prevelance of the phenotype (if binary, otherwise NA)
- mean: The phenotype mean in the general population (if continuous, otherwise NA)
- sd: The phenotype sd in the general population (if continuous, otherwise NA)
- label: A human readable name for the GWAS phenotype. Wrap in quotes if multiple words. For example, "Body Mass Index".

#### GWAS summary statistic format

The following column names are expected in the GWAS summary statistics files:

- SNP: RSID for variant (required)
- A1: Allele 1 (effect allele) (required)
- A2: Allele 2 (required)
- P: P-value of association (required)
- OR: Odds ratio effect size (required if binary)
- BETA: BETA effect size (required if continuous)
- SE: Standard error of log(OR) or BETA (required)
- N: Total sample size (required)
- FREQ: Allele frequency in GWAS sample (optional)
- INFO: Imputation quality (optional)


***

Please contact Oliver Pain (oliver.pain@kcl.ac.uk) for any questions or comments.

***
