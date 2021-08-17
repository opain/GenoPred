# Genotype-based Prediction Pipeline (GenoPredPipe)

This is a snakemake pipeline for running the GenoPred scripts. Currently it can intake two types of data, 23andMe data for an individual, or PLINK binary format (.bed/.bim/.fam) for any sample. The pipeline currently identifies the ancestry of each individual, calculates ancestry-matched reference-projected principal components, and calculates polygenic scores using a range of methods standardised using an ancestry matched reference. Finally, the pipeline provides a report of the results either for each individual or the sample overall.

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

### Step 1

Download test data.

```bash
wget -O test_data.tar.gz https://zenodo.org/record/5205711/files/test_data.tar.gz?download=1
tar -xf test_data.tar.gz
rm test_data.tar.gz
```

### Step 2

Run full pipeline.

```bash
snakemake --profile slurm --use-conda run_create_reports
```



Note. alternative GWAS and target sample datasets can be specified by modifying the config.yaml file.

Please contact Oliver Pain (oliver.pain@kcl.ac.uk) for any questions or comments.