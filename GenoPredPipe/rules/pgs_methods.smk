# Create PC score files specific to each population
rule pop_pc_scoring:
  input:
    rules.get_dependencies.output,
    "../Scripts/ancestry_score_file_creator/ancestry_score_file_creator.R",
    "../Scripts/functions/misc.R"
  output:
    "resources/data/ref/pc_score_files/{population}/ref.{population}.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/ancestry_score_file_creator/ancestry_score_file_creator.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{wildcards.population}.keep \
      --plink plink \
      --plink2 plink2 \
      --n_pcs 6 \
      --output resources/data/ref/pc_score_files/{wildcards.population}/ref.{wildcards.population}"

populations=["AFR","AMR","EAS","EUR","SAS"]

rule run_pop_pc_scoring:
  input: expand("resources/data/ref/pc_score_files/{population}/ref.{population}.scale", population=populations)
  
##
# QC and format GWAS summary statistics
##

import pandas as pd
gwas_list_df = pd.read_table(config["gwas_list"], sep=r'\s+')
gwas_list_df_eur = gwas_list_df.loc[gwas_list_df['population'] == 'EUR']

score_list_file = Path(config["score_list"])
if score_list_file.is_file():
  score_list_df = pd.read_table(config["score_list"], sep=r'\s+')
else:
  score_list_df = pd.DataFrame(columns = ["name", "path", "population", "sampling", "prevalence", "mean", "sd", "label"])

rule sumstat_prep:
  input:
    config['gwas_list'],
    rules.get_dependencies.output,
    "../Scripts/sumstat_cleaner/sumstat_cleaner.R"
  output:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    path= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0]
  shell:
    "Rscript ../Scripts/sumstat_cleaner/sumstat_cleaner.R \
      --sumstats {params.path} \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_freq_chr resources/data/ref/freq_files/{params.population}/ref.{params.population}.chr \
      --output resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned"
    
rule run_sumstat_prep:
  input: expand("resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz", gwas=gwas_list_df['name'])

##
# pT+clump (sparse, nested)
##

rule prs_scoring_ptclump:
  input:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    "../Scripts/pgs_methods/ptclump.R",
    "../Scripts/functions/misc.R"
  output:
    "resources/data/ref/prs_score_files/ptclump/{gwas}/ref.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/ptclump.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{params.population}.keep \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink1 plink \
      --plink2 plink2 \
      --output resources/data/ref/prs_score_files/ptclump/{wildcards.gwas}/ref.{wildcards.gwas} \
      --ref_pop_scale resources/data/ref/ref.keep.list \
      --test {params.testing}"
  
rule run_prs_scoring_ptclump:
  input: expand("resources/data/ref/prs_score_files/ptclump/{gwas}/ref.{gwas}.EUR.scale", gwas=gwas_list_df['name'])

##
# DBSLMM
##

rule prs_scoring_dbslmm:
  input:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.get_dependencies.output,
    "../Scripts/pgs_methods/dbslmm.R",
    "../Scripts/functions/misc.R"
  output:
    "resources/data/ref/prs_score_files/dbslmm/{gwas}/ref.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    sampling= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'sampling'].iloc[0],
    prevalence= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'prevalence'].iloc[0],
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/dbslmm.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{params.population}.keep \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink resources/software/plink/plink \
      --ld_blocks resources/data/ld_blocks/{params.population} \
      --rscript Rscript \
      --dbslmm resources/software/dbslmm/software \
      --munge_sumstats resources/software/ldsc/munge_sumstats.py \
      --ldsc resources/software/ldsc/ldsc.py \
      --ldsc_ref resources/data/ldsc_ref/eur_w_ld_chr \
      --hm3_snplist resources/data/hm3_snplist/w_hm3.snplist \
      --sample_prev {params.sampling} \
      --pop_prev {params.prevalence} \
      --output resources/data/ref/prs_score_files/dbslmm/{wildcards.gwas}/ref.{wildcards.gwas} \
      --ref_pop_scale resources/data/ref/ref.keep.list \
      --test {params.testing}"

rule run_prs_scoring_dbslmm:
  input: expand("resources/data/ref/prs_score_files/dbslmm/{gwas}/ref.{gwas}.EUR.scale", gwas=gwas_list_df['name'])

##
# PRScs
##

rule prs_scoring_prscs:
  resources: 
    mem_mb=80000,
    cpus=10,
    time_min=800
  input:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.get_dependencies.output,
    "../Scripts/pgs_methods/prscs.R",
    "../Scripts/functions/misc.R",
    rules.download_prscs_software.output,
    rules.download_prscs_ref_1kg_eur.output
  output:
    "resources/data/ref/prs_score_files/prscs/{gwas}/ref.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    testing=config["testing"]
  shell:
    "export MKL_NUM_THREADS=1; \
     export NUMEXPR_NUM_THREADS=1; \
     export OMP_NUM_THREADS=1; \
     Rscript ../Scripts/pgs_methods/prscs.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink2 plink2 \
      --memory 5000 \
      --output resources/data/ref/prs_score_files/prscs/{wildcards.gwas}/ref.{wildcards.gwas} \
      --ref_pop_scale resources/data/ref/ref.keep.list \
      --PRScs_path resources/software/prscs/PRScs.py \
      --PRScs_ref_path resources/data/prscs_ref/ldblk_1kg_eur \
      --n_cores {resources.cpus} \
      --phi_param 1e-6,1e-4,1e-2,1,auto \
      --seed 1 \
      --test {params.testing}"

rule run_prs_scoring_prscs:
  input: expand("resources/data/ref/prs_score_files/prscs/{gwas}/ref.{gwas}.EUR.scale", gwas=gwas_list_df_eur['name'])

##
# SBayesR
##

rule prs_scoring_sbayesr:
  resources: 
    mem_mb=80000,
    cpus=10
  input:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.get_dependencies.output,
    "../Scripts/pgs_methods/sbayesr.R",
    "../Scripts/functions/misc.R",
    rules.download_gctb_ref.output,
    rules.download_gctb_software.output
  output:
    "resources/data/ref/prs_score_files/sbayesr/{gwas}/ref.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/sbayesr.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink plink \
      --gctb resources/software/gctb/gctb_2.03beta_Linux/gctb \
      --ld_matrix_chr resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr \
      --memory {resources.mem_mb} \
      --robust T \
      --n_cores {resources.cpus} \
      --output resources/data/ref/prs_score_files/sbayesr/{wildcards.gwas}/ref.{wildcards.gwas} \
      --ref_pop_scale resources/data/ref/ref.keep.list \
      --test {params.testing}"

rule run_prs_scoring_sbayesr:
  input: expand("resources/data/ref/prs_score_files/sbayesr/{gwas}/ref.{gwas}.EUR.scale", gwas=gwas_list_df_eur['name'])

##
# lassosum
##

rule prs_scoring_lassosum:
  input:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.get_dependencies.output,
    "../Scripts/pgs_methods/lassosum.R",
    "../Scripts/functions/misc.R"
  output:
    "resources/data/ref/prs_score_files/lassosum/{gwas}/ref.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/lassosum.R \
     --ref_plink_chr resources/data/ref/ref.chr \
     --ref_keep resources/data/ref/keep_files/{params.population}.keep \
     --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
     --output resources/data/ref/prs_score_files/lassosum/{wildcards.gwas}/ref.{wildcards.gwas} \
     --plink2 plink2 \
     --ref_pop_scale resources/data/ref/ref.keep.list \
      --test {params.testing}"
    
rule run_prs_scoring_lassosum:
  input: expand("resources/data/ref/prs_score_files/lassosum/{gwas}/ref.{gwas}.EUR.scale", gwas=gwas_list_df['name'])

##
# LDpred2
##

rule prs_scoring_ldpred2:
  resources: 
    mem_mb=80000,
    cpus=10,
    time_min=800
  input:
    rules.get_dependencies.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    "../Scripts/pgs_methods/ldpred2.R",
    "../Scripts/functions/misc.R",
    rules.download_ldpred2_ref.output
  output:
    "resources/data/ref/prs_score_files/ldpred2/{gwas}/ref.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/ldpred2.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/EUR.keep \
      --ldpred2_ref_dir resources/data/ldpred2_ref \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink2 plink2 \
      --memory {resources.mem_mb} \
      --n_cores {resources.cpus} \
      --output resources/data/ref/prs_score_files/ldpred2/{wildcards.gwas}/ref.{wildcards.gwas} \
      --ref_pop_scale resources/data/ref/ref.keep.list \
      --test {params.testing}"
    
rule run_prs_scoring_ldpred2:
  input: expand("resources/data/ref/prs_score_files/ldpred2/{gwas}/ref.{gwas}.EUR.scale", gwas=gwas_list_df_eur['name'])

##
# LDAK MegaPRS
##

rule prs_scoring_megaprs:
  resources: 
    mem_mb=20000,
    cpus=5,
    time_min=800
  input:
    rules.get_dependencies.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    "../Scripts/pgs_methods/megaprs.R",
    "../Scripts/functions/misc.R",
    rules.download_ldak_highld.output,
    rules.download_ldak.output,
    rules.download_ldak_bld.output
  output:
    "resources/data/ref/prs_score_files/megaprs/{gwas}/ref.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/megaprs.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{params.population}.keep \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink1 plink \
      --plink2 plink2 \
      --ldak resources/software/ldak/ldak5.1.linux \
      --ldak_map resources/data/ldak_map/genetic_map_b37 \
      --ldak_tag resources/data/ldak_bld \
      --ldak_highld resources/data/ldak_highld/highld.txt \
      --memory {resources.mem_mb} \
      --n_cores {resources.cpus} \
      --output resources/data/ref/prs_score_files/megaprs/{wildcards.gwas}/ref.{wildcards.gwas} \
      --ref_pop_scale resources/data/ref/ref.keep.list \
      --test {params.testing}"
    
rule run_prs_scoring_megaprs:
  input: expand("resources/data/ref/prs_score_files/megaprs/{gwas}/ref.{gwas}.EUR.scale", gwas=gwas_list_df['name'])

##
# Process externally created score files
##

# Read in list of external score files
score_list_file = Path(config["score_list"])
if score_list_file.is_file():
  score_list_df = pd.read_table(config["score_list"], sep=r'\s+')
else:
  score_list_df = pd.DataFrame(columns = ["name", "path", "population", "sampling", "prevalence", "mean", "sd", "label"])

rule prs_scoring_external:
  input:
    config['score_list'],
    rules.get_dependencies.output,
    lambda w: score_list_df.loc[score_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0],
    "../Scripts/external_score_processor/external_score_processor.R",
    "../Scripts/functions/misc.R"
  output:
    "resources/data/ref/prs_score_files/external/{gwas}/ref.{gwas}.EUR.scale"
  params:
    score= lambda w: score_list_df.loc[score_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0],
    population= lambda w: score_list_df.loc[score_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0]
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/external_score_processor/external_score_processor.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --score {params.score} \
      --plink2 plink2 \
      --output resources/data/ref/prs_score_files/external/{wildcards.gwas}/ref.{wildcards.gwas} \
      --ref_pop_scale resources/data/ref/ref.keep.list"
    
rule run_prs_scoring_external:
  input: expand("resources/data/ref/prs_score_files/external/{gwas}/ref.{gwas}.EUR.scale", gwas=score_list_df['name'])

##
# Estimate R2/AUC of PRS using lassosum pseudovalidate
##

rule pseudovalidate_prs:
  input:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.get_dependencies.output,
    "../Scripts/lassosum_pseudovalidate/lassosum_pseudovalidate.R"
  output:
    "resources/data/ref/prs_pseudoval/{gwas}/lassosum_pseudo_{gwas}.pseudovalidate.png"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
  shell:
    "Rscript ../Scripts/lassosum_pseudovalidate/lassosum_pseudovalidate.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{params.population}.keep \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --prune_mhc T \
      --output resources/data/ref/prs_pseudoval/{wildcards.gwas}/lassosum_pseudo_{wildcards.gwas} \
      --plink plink \
      --n_cores 1"

rule run_pseudovalidate_prs:
  input: expand("resources/data/ref/prs_pseudoval/{gwas}/lassosum_pseudo_{gwas}.pseudovalidate.png", gwas=gwas_list_df['name'])

