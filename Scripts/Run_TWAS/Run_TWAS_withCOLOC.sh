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
    --fusion_dir)
    fusion_dir="$2"
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

module add apps/R/3.6.0

gwas=$(cut -d' ' -f 2 $gwas_list | awk "NR==${SLURM_ARRAY_TASK_ID}")
name=$(cut -d' ' -f 1 $gwas_list | awk "NR==${SLURM_ARRAY_TASK_ID}")
out=$(echo ${outdir}/${name}_withCOLOC/${name})

Rscript /users/k1806347/brc_scratch/Software/MyGit/GenoPred/Scripts/Run_TWAS/Run_TWAS.R \
  --gwas ${gwas} \
  --pos ${pos} \
  --coloc_P 0.05 \
  --out ${out} \
  --fusion_dir ${fusion_dir} \
  --ncores ${ncores}
