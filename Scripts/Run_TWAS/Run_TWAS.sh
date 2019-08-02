#!/bin/bash
#$ -b y

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    --gwas_list)
    gwas_list="$2"
    shift # past argument
    shift # past value
    ;;
    --pos)
    pos="$2"
    shift # past argument
    shift # past value
    ;;
    --outdir)
    outdir="$2"
    shift # past argument
    shift # past value
    ;;
    --ncores)
    ncores="$2"
    shift # past argument
    shift # past value
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

module add general/R/3.5.0

gwas=$(cut -d' ' -f 2 $gwas_list | awk "NR==${SGE_TASK_ID}")
name=$(cut -d' ' -f 1 $gwas_list | awk "NR==${SGE_TASK_ID}")
out=$(echo ${outdir}/${name}/${name})

Rscript /users/k1806347/brc_scratch/Software/Scripts/GenoPred_Pipeline/Run_TWAS.R \
  --gwas ${gwas} \
  --pos ${pos} \
  --out ${out} \
  --ncores ${ncores}
