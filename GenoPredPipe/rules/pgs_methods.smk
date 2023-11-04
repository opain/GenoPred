# Create PC score files specific to each super population
rule super_pop_pc_scoring:
  input:
    ref=rules.prep_1kg.output
  output:
    "resources/data/1kg/pc_score_files/{population}/1KGPhase3.w_hm3.{population}.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/ancestry_score_file_creator/ancestry_score_file_creator.R \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --ref_keep resources/data/1kg/keep_files/{wildcards.population}_samples.keep \
      --plink plink \
      --plink2 plink2 \
      --n_pcs 100 \
      --output resources/data/1kg/pc_score_files/{wildcards.population}/1KGPhase3.w_hm3.{wildcards.population}"

populations=["AFR","AMR","EAS","EUR","SAS"]

rule run_super_pop_pc_scoring:
  input: expand("resources/data/1kg/pc_score_files/{population}/1KGPhase3.w_hm3.{population}.scale", population=populations)
  
##
# QC and format GWAS summary statistics
##

import pandas as pd
gwas_list_df = pd.read_table(config["gwas_list"], sep=' ')
gwas_list_df_eur = gwas_list_df.loc[gwas_list_df['population'] == 'EUR']

rule sumstat_prep:
  input:
    config['gwas_list'],
    rules.prep_1kg.output
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
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --ref_freq_chr resources/data/1kg/freq_files/{params.population}/1KGPhase3.w_hm3.{params.population}.chr \
      --output resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned"
    
rule run_sumstat_prep:
  input: expand("resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz", gwas=gwas_list_df['name'])

##
# pT+clump (sparse, nested)
##

rule prs_scoring_pt_clump:
  input:
    rules.prep_1kg.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz"
  output:
    "resources/data/1kg/prs_score_files/pt_clump/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
  shell:
    "Rscript ../Scripts/polygenic_score_file_creator/polygenic_score_file_creator_plink2.R \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink1 plink \
      --plink2 plink2 \
      --output resources/data/1kg/prs_score_files/pt_clump/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list"
  
rule run_prs_scoring_pt_clump:
  input: expand("resources/data/1kg/prs_score_files/pt_clump/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df['name'])

##
# DBSLMM
##

rule prs_scoring_dbslmm:
  input:
    rules.prep_1kg.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.download_ldsc.output,
    rules.dowload_ldsc_ref.output,
    rules.download_dbslmm.output,
    rules.download_ld_blocks.output,
    rules.download_plink.output
  output:
    "resources/data/1kg/prs_score_files/dbslmm/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    sampling= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'sampling'].iloc[0],
    prevalence= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'prevalence'].iloc[0],
  shell:
    "Rscript ../Scripts/polygenic_score_file_creator_DBSLMM/polygenic_score_file_creator_DBSLMM_plink2.R \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep \
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
      --output resources/data/1kg/prs_score_files/dbslmm/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list"

rule run_prs_scoring_dbslmm:
  input: expand("resources/data/1kg/prs_score_files/dbslmm/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df['name'])

##
# PRScs
##

rule prs_scoring_prscs:
  resources: 
    mem_mb=80000,
    cpus=10,
    time_min=800
  input:
    rules.prep_1kg.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.download_plink.output,
    rules.download_prscs_ref_1kg_eur.output,
    rules.download_prscs_software.output
  output:
    "resources/data/1kg/prs_score_files/prscs/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "export MKL_NUM_THREADS=1; \
     export NUMEXPR_NUM_THREADS=1; \
     export OMP_NUM_THREADS=1; \
     Rscript ../Scripts/polygenic_score_file_creator_PRScs/polygenic_score_file_creator_PRScs_plink2.R \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink2 plink2 \
      --memory 5000 \
      --output resources/data/1kg/prs_score_files/prscs/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list \
      --PRScs_path resources/software/prscs/PRScs.py \
      --PRScs_ref_path resources/data/prscs_ref/ldblk_1kg_eur \
      --n_cores {resources.cpus} \
      --phi_param 1e-6,1e-4,1e-2,1,auto \
      --seed 1"

rule run_prs_scoring_prscs:
  input: expand("resources/data/1kg/prs_score_files/prscs/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df_eur['name'])

##
# SBayesR
##

rule prs_scoring_sbayesr:
  resources: 
    mem_mb=80000,
    cpus=10
  input:
    rules.prep_1kg.output,
    rules.merge_1kg_GW.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.download_plink.output,
    rules.download_gctb_ref.output,
    rules.download_gctb_software.output
  output:
    "resources/data/1kg/prs_score_files/sbayesr/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/polygenic_score_file_creator_SBayesR/polygenic_score_file_creator_SBayesR_plink2.R \
      --ref_plink resources/data/1kg/1KGPhase3.w_hm3.GW \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink plink \
      --gctb resources/software/gctb/gctb_2.03beta_Linux/gctb \
      --ld_matrix_chr resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr \
      --memory {resources.mem_mb} \
      --robust T \
      --n_cores {resources.cpus} \
      --output resources/data/1kg/prs_score_files/sbayesr/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list"

rule run_prs_scoring_sbayesr:
  input: expand("resources/data/1kg/prs_score_files/sbayesr/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df_eur['name'])

##
# lassosum
##

rule prs_scoring_lassosum:
  input:
    rules.prep_1kg.output,
    rules.merge_1kg_GW.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.download_plink.output,
    rules.install_lassosum.output
  output:
    "resources/data/1kg/prs_score_files/lassosum/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
  shell:
    "Rscript ../Scripts/polygenic_score_file_creator_lassosum/polygenic_score_file_creator_lassosum_plink2.R \
     --ref_plink_gw resources/data/1kg/1KGPhase3.w_hm3.GW \
     --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep \
     --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
     --output resources/data/1kg/prs_score_files/lassosum/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
     --plink2 plink2 \
     --ref_pop_scale resources/data/1kg/super_pop_keep.list"
    
rule run_prs_scoring_lassosum:
  input: expand("resources/data/1kg/prs_score_files/lassosum/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df['name'])

##
# LDpred2
##

rule prs_scoring_ldpred2:
  resources: 
    mem_mb=80000,
    cpus=10,
    time_min=800
  input:
    rules.prep_1kg.output,
    rules.merge_1kg_GW.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.download_ldpred2_ref.output
  output:
    "resources/data/1kg/prs_score_files/ldpred2/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/polygenic_score_file_creator_LDPred2/polygenic_score_file_creator_LDPred2_LDPredRef_plink2.R \
      --ref_plink resources/data/1kg/1KGPhase3.w_hm3.GW \
      --ref_keep resources/data/1kg/keep_files/EUR_samples.keep \
      --ldpred2_ref_dir resources/data/ldpred2_ref \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink2 plink2 \
      --memory {resources.mem_mb} \
      --n_cores {resources.cpus} \
      --output resources/data/1kg/prs_score_files/ldpred2/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list"
    
rule run_prs_scoring_ldpred2:
  input: expand("resources/data/1kg/prs_score_files/ldpred2/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df_eur['name'])

##
# LDAK MegaPRS
##

rule prs_scoring_megaprs:
  resources: 
    mem_mb=20000,
    cpus=5,
    time_min=800
  input:
    rules.prep_1kg.output,
    rules.merge_1kg_GW.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.download_ldak.output,
    rules.download_ldak_map.output,
    rules.download_ldak_bld.output,
    rules.download_ldak_highld.output
  output:
    "resources/data/1kg/prs_score_files/megaprs/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
  shell:
    "Rscript ../Scripts/ldak_mega_prs/ldak_mega_prs.R \
      --ref_plink resources/data/1kg/1KGPhase3.w_hm3.GW \
      --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --plink1 plink \
      --plink2 plink2 \
      --ldak resources/software/ldak/ldak5.1.linux \
      --ldak_map resources/data/ldak_map/genetic_map_b37 \
      --ldak_tag resources/data/ldak_bld \
      --ldak_highld resources/data/ldak_highld/highld.txt \
      --memory {resources.mem_mb} \
      --n_cores {resources.cpus} \
      --output resources/data/1kg/prs_score_files/megaprs/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list"
    
rule run_prs_scoring_megaprs:
  input: expand("resources/data/1kg/prs_score_files/megaprs/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=gwas_list_df['name'])

##
# Process externally created score files
##

# Read in list of external score files
score_list_file = Path(config["score_list"])
if score_list_file.is_file():
  score_list_df = pd.read_table(config["score_list"], sep=' ')
else:
  score_list_df = pd.DataFrame(columns = ["name", "path", "population", "sampling", "prevalence", "mean", "sd", "label"])

rule prs_scoring_external:
  input:
    config['score_list'],
    rules.prep_1kg.output,
    lambda w: score_list_df.loc[score_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0],
  output:
    "resources/data/1kg/prs_score_files/external/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale"
  params:
    score= lambda w: score_list_df.loc[score_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0],
    population= lambda w: score_list_df.loc[score_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0]
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/external_score_processor/external_score_processor_plink2.R \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --score {params.score} \
      --plink2 plink2 \
      --output resources/data/1kg/prs_score_files/external/{wildcards.gwas}/1KGPhase3.w_hm3.{wildcards.gwas} \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list"
    
rule run_prs_scoring_external:
  input: expand("resources/data/1kg/prs_score_files/external/{gwas}/1KGPhase3.w_hm3.{gwas}.EUR.scale", gwas=score_list_df['name'])

##
# Estimate R2/AUC of PRS using lassosum pseudovalidate
##

rule pseudovalidate_prs:
  input:
    rules.merge_1kg_GW.output,
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.install_lassosum.output
  output:
    "resources/data/1kg/prs_pseudoval/{gwas}/lassosum_pseudo_{gwas}.pseudovalidate.png"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
  shell:
    "Rscript ../Scripts/lassosum_pseudovalidate/lassosum_pseudovalidate.R \
      --ref_plink_gw resources/data/1kg/1KGPhase3.w_hm3.GW \
      --ref_keep resources/data/1kg/keep_files/{params.population}_samples.keep \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --prune_mhc T \
      --output resources/data/1kg/prs_pseudoval/{wildcards.gwas}/lassosum_pseudo_{wildcards.gwas} \
      --plink plink \
      --n_cores 1"

rule run_pseudovalidate_prs:
  input: expand("resources/data/1kg/prs_pseudoval/{gwas}/lassosum_pseudo_{gwas}.pseudovalidate.png", gwas=gwas_list_df['name'])

rule pipeline_prep:
  input:
    rules.run_super_pop_pc_scoring.input,
    rules.run_prs_scoring_pt_clump.input,
    rules.run_prs_scoring_dbslmm.input,
    rules.run_pseudovalidate_prs.input
