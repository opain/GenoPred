#!/bin/bash

# Find conda.sh
CONDA_EXE=$(which conda)
CONDA_BASE=$(dirname $(dirname $CONDA_EXE))

# Source the conda.sh script
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate genopred

# Run Snakemake
snakemake $*
