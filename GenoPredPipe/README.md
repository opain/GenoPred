# GenoPredPipe

This is a snakemake pipeline for running the GenoPred scripts.

## Getting started

### Step 1

Clone the repository

```
git clone git@github.com:https://github.com/opain/GenoPred.git
cd GenoPred/GenoPredPipe
```

### Step 2

Install [Anaconda](https://conda.io/en/latest/miniconda.html).

**Linux:**
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
```

Install Python 3.8, Snakemake 5.32, and the basic project dependencies

```
conda activate base
conda install python=3.8
conda install -c conda-forge mamba
mamba install -c bioconda -c conda-forge snakemake-minimal==5.32.2
```

## Running pipeline

### Download test data

```
wget -O test_data.tar.gz --no-check-certificate 'https://docs.google.com/uc?export=download&id=1Gqto9__A8eH-SrSWvVjT_XexXFlsKJMu'
tar -xf test_data.tar.gz
```
### Run full pipeline

```
snakemake --profile slurm --use-conda run_create_reports
```


