import pandas as pd
target_list_df = pd.read_table(config["target_list"], sep=' ')

target_list_df_23andMe = target_list_df.loc[target_list_df['type'] == '23andMe']
target_list_df_samp_imp_plink1 = target_list_df.loc[target_list_df['type'] == 'samp_imp_plink1']
target_list_df_samp_imp_bgen = target_list_df.loc[target_list_df['type'] == 'samp_imp_bgen']
target_list_df_samp_imp_vcf = target_list_df.loc[target_list_df['type'] == 'samp_imp_vcf']
target_list_df_samp_imp = target_list_df.loc[(target_list_df['type'] == 'samp_imp_plink1') | (target_list_df['type'] == 'samp_imp_bgen') | (target_list_df['type'] == 'samp_imp_vcf')]
target_list_df_samp_imp_indiv_report = target_list_df_samp_imp.loc[(target_list_df['indiv_report'] == 'T')]

####
# Format target data
####

##
# 23andMe
##
# Largely based on Impute.Me by Lasse Folkersen

rule format_impute_23andme_target:
  resources: 
    mem_mb=90000,
    cpus=6
  input:
    config['target_list'],
    rules.download_impute2_data.output,
    rules.download_qctool2.output,
    "../Scripts/23andMe_imputer/23andMe_imputer.R"
  output:
    touch("resources/data/target_checks/{name}/format_impute_23andme_target.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    name= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'name'].iloc[0],
    path= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0],
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0]
  shell:
    "Rscript ../Scripts/23andMe_imputer/23andMe_imputer.R \
      --FID {params.name} \
      --IID {params.name} \
      --geno {params.path} \
      --plink plink \
      --out {params.output}/{params.name}/{params.name} \
      --ref resources/data/impute2/1000GP_Phase3 \
      --shapeit shapeit \
      --impute2 impute2 \
      --n_core {resources.cpus} \
      --qctool resources/software/qctool2/qctool"

rule run_format_impute_23andme_target:
  input: expand("resources/data/target_checks/{name}/format_impute_23andme_target_{name}.done", name=target_list_df_23andMe['name'])

##
# Formatting without imputation
##
# Here we will use a script that can accept multiple genetic data formats, insert RSIDs if not present, extract HapMap3 SNPs and output plink hard-call binaries.

rule format_target:
  input:
    config['target_list'],
    rules.prep_1kg.output,
    rules.install_liftover.output,
    rules.download_liftover_track.output,
    "../Scripts/geno_to_plink/geno_to_plink.R"
  output:
    touch("resources/data/target_checks/{name}/format_target_{chr}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    path= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0],
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0],
    type= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'type'].iloc[0]
  shell:
    "Rscript ../Scripts/geno_to_plink/geno_to_plink.R \
      --target {params.path}.chr{wildcards.chr} \
      --format {params.type} \
      --ref resources/data/1kg/1KGPhase3.w_hm3.chr{wildcards.chr} \
      --plink2 plink2 \
      --liftover resources/software/liftover/liftover \
      --liftover_track resources/data/liftover/hg19ToHg38.over.chain.gz \
      --out {params.output}/{wildcards.name}/{wildcards.name}.hm3.chr{wildcards.chr}"

rule run_format_target:
  input: 
    lambda w: expand("resources/data/target_checks/{name}/format_target_{chr}.done", name=w.name, chr=range(1, 23))
  output:
    touch("resources/data/target_checks/{name}/format_target.done")

rule run_format_target_2:
  input: 
    expand("resources/data/target_checks/{name}/format_target.done", name=target_list_df_samp_imp['name'])

####
# Harmonise with reference
####

rule harmonise_target_with_ref:
  input:
    rules.prep_1kg.output,
    lambda w: "resources/data/target_checks/" + w.name + "/format_impute_23andme_target.done" if target_list_df.loc[target_list_df['name'] == w.name, 'type'].iloc[0] == '23andMe' else ("resources/data/target_checks/{name}/format_target.done" if target_list_df.loc[target_list_df['name'] == w.name, 'type'].iloc[0] == 'samp_imp_plink1' else ("resources/data/target_checks/{name}/format_target.done" if target_list_df.loc[target_list_df['name'] == w.name, 'type'].iloc[0] == 'samp_imp_bgen' else "resources/data/target_checks/{name}/format_target.done")),
    "../Scripts/hm3_harmoniser/hm3_harmoniser.R"
  output:
    touch("resources/data/target_checks/{name}/harmonise_target_with_ref.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0]
  shell:
    "Rscript ../Scripts/hm3_harmoniser/hm3_harmoniser.R \
      --target {params.output}/{wildcards.name}/{wildcards.name}.hm3.chr \
      --ref resources/data/1kg/1KGPhase3.w_hm3.chr \
      --plink plink \
      --out {params.output}/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr"

rule run_harmonise_target_with_ref:
  input: expand("resources/data/target_checks/{name}/harmonise_target_with_ref.done", name=target_list_df['name'])
  
# Delete temporary files
rule delete_temp_target_samp_files:
  input:
    "resources/data/target_checks/{name}/harmonise_target_with_ref.done"
  output:
    touch("resources/data/target_checks/{name}/delete_temp_target_samp_files.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0]
  shell:
    "rm {params.output}/{wildcards.name}/{wildcards.name}.hm3.chr*"

####
# Identify super_population
####

rule target_super_population:
  input:
    "resources/data/target_checks/{name}/harmonise_target_with_ref.done",
    lambda w: "resources/data/target_checks/" + w.name + "/harmonise_target_with_ref.done" if target_list_df.loc[target_list_df['name'] == w.name, 'type'].iloc[0] == '23andMe' else "resources/data/target_checks/" + w.name + "/delete_temp_target_samp_files.done",
    "../Scripts/Ancestry_identifier/Ancestry_identifier.R"
  output:
    touch("resources/data/target_checks/{name}/target_super_population.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0]
  shell:
    "Rscript ../Scripts/Ancestry_identifier/Ancestry_identifier.R \
      --target_plink_chr {params.output}/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr \
      --ref_plink_chr resources/data/1kg/1KGPhase3.w_hm3.chr \
      --n_pcs 6 \
      --plink plink \
      --plink2 plink2 \
      --output {params.output}/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry \
      --ref_pop_scale resources/data/1kg/super_pop_keep.list \
      --pop_data resources/data/1kg/ref_pop_dat.txt \
      --prob_thresh 0.5"

rule run_target_super_population:
  input: expand("resources/data/target_checks/{name}/target_super_population.done", name=target_list_df['name'])

# Create a file listing target samples and super population assignments
checkpoint ancestry_reporter:
  input:
    "resources/data/target_checks/{name}/target_super_population.done",
    "scripts/ancestry_reporter.R"
  output:
    touch("resources/data/target_checks/{name}/ancestry_reporter.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0]
  shell:
    "Rscript scripts/ancestry_reporter.R {wildcards.name} {params.output}"

rule run_ancestry_reporter:
  input: expand("resources/data/target_checks/{name}/ancestry_reporter.done", name=target_list_df['name'])

####
# Super population outlier detection and target sample specific PC calculation
####

rule target_super_population_outlier_detection:
  resources: 
    mem_mb=15000
  input:
    "resources/data/target_checks/{name}/ancestry_reporter.done",
    "../Scripts/Population_outlier/Population_outlier.R"
  output:
    touch('resources/data/target_checks/{name}/target_super_population_outlier_detection.done')
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0]
  shell:
    "ls {params.output}/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.*.keep > {params.output}/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.keep_list; Rscript ../Scripts/Population_outlier/Population_outlier.R \
      --target_plink_chr {params.output}/{wildcards.name}/{wildcards.name}.1KGphase3.hm3.chr \
      --target_keep {params.output}/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.keep_list \
      --n_pcs 10 \
      --maf 0.05 \
      --geno 0.02 \
      --hwe 1e-6 \
      --memory {resources.mem_mb} \
      --plink plink \
      --plink2 plink2 \
      --output {params.output}/{wildcards.name}/ancestry/outlier_detection/{wildcards.name}.outlier_detection"

# Create a function summarising which populations are present in target      
from pathlib import Path

def ancestry_munge(x):
    checkpoint_output = checkpoints.ancestry_reporter.get(name=x).output[0]
    checkpoint_output = target_list_df.loc[target_list_df['name'] == x, 'output'].iloc[0] + "/" + x + "/ancestry/ancestry_report.txt"
    ancestry_report_df = pd.read_table(checkpoint_output, sep=' ')
    return ancestry_report_df['population'].tolist()
    