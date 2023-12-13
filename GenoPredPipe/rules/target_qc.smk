import pandas as pd
target_list_df = pd.read_table(config["target_list"], sep=r'\s+')

target_list_df_23andMe = target_list_df.loc[target_list_df['type'] == '23andMe']
samp_types = ['samp_imp_plink1', 'samp_imp_bgen', 'samp_imp_vcf']
target_list_df_samp_imp = target_list_df[target_list_df['type'].isin(samp_types)]
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
    touch("{outdir}/resources/data/target_checks/{name}/format_impute_23andme_target.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    name= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'name'].iloc[0],
    path= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0]
  shell:
    "Rscript ../Scripts/23andMe_imputer/23andMe_imputer.R \
      --FID {params.name} \
      --IID {params.name} \
      --geno {params.path} \
      --plink plink \
      --out {outdir}/{params.name}/geno/{params.name} \
      --ref resources/data/impute2/1000GP_Phase3 \
      --shapeit shapeit \
      --impute2 impute2 \
      --n_core {resources.cpus} \
      --qctool resources/software/qctool2/qctool"

rule run_format_impute_23andme_target:
  input: expand("{outdir}/resources/data/target_checks/{name}/format_impute_23andme_target_{name}.done", name=target_list_df_23andMe['name'], outdir=outdir)

##
# Convert to PLINK and harmonise with reference
##

rule format_target:
  input:
    config['target_list'],
    rules.get_dependencies.output,
    "../Scripts/format_target/format_target.R"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/format_target_{chr}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    path= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0],
    type= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'type'].iloc[0]
  shell:
    "Rscript ../Scripts/format_target/format_target.R \
      --target {params.path}.chr{wildcards.chr} \
      --format {params.type} \
      --ref resources/data/ref/ref.chr{wildcards.chr} \
      --output {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr{wildcards.chr}"

rule run_format_target:
  input: 
    lambda w: expand("{outdir}/resources/data/target_checks/{name}/format_target_{chr}.done", name=w.name, chr=range(1, 23), outdir=outdir)
  output:
    touch("{outdir}/resources/data/target_checks/{name}/format_target.done")

####
# Population classification
####

rule ancestry_inference:
  input:
    "{outdir}/resources/data/target_checks/{name}/format_target.done",
    "../Scripts/Ancestry_identifier/Ancestry_identifier.R"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/ancestry_inference.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/Ancestry_identifier/Ancestry_identifier.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --ref_plink_chr resources/data/ref/ref.chr \
      --output {outdir}/{wildcards.name}/ancestry/{wildcards.name}.Ancestry \
      --pop_data resources/data/ref/ref.pop.txt"

rule run_ancestry_inference:
  input: expand("{outdir}/resources/data/target_checks/{name}/ancestry_inference.done", name=target_list_df['name'], outdir=outdir)

# Create a file listing target samples and population assignments
checkpoint ancestry_reporter:
  input:
    "{outdir}/resources/data/target_checks/{name}/ancestry_inference.done",
    "scripts/ancestry_reporter.R"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/ancestry_reporter.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript scripts/ancestry_reporter.R {wildcards.name} {outdir}"

rule run_ancestry_reporter:
  input: expand("{outdir}/resources/data/target_checks/{name}/ancestry_reporter.done", name=target_list_df['name'], outdir=outdir)

####
# Population outlier detection and target sample specific PC calculation
####

rule outlier_detection:
  resources: 
    mem_mb=15000
  input:
    "{outdir}/resources/data/target_checks/{name}/ancestry_reporter.done",
    "../Scripts/Population_outlier/Population_outlier.R"
  output:
    touch('{outdir}/resources/data/target_checks/{name}/outlier_detection.done')
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "ls {outdir}/{wildcards.name}/ancestry/{wildcards.name}.Ancestry.model_pred.*.keep > {outdir}/{wildcards.name}/ancestry/{wildcards.name}.Ancestry.model_pred.keep_list; Rscript ../Scripts/Population_outlier/Population_outlier.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --target_keep {outdir}/{wildcards.name}/ancestry/keep_files/model_based/keep_list.txt \
      --n_pcs 10 \
      --maf 0.05 \
      --geno 0.02 \
      --hwe 1e-6 \
      --memory {resources.mem_mb} \
      --plink plink \
      --plink2 plink2 \
      --output {outdir}/{wildcards.name}/pcs/within_sample/{wildcards.name}.outlier_detection"

# Create a function summarising which populations are present in target      
from pathlib import Path

def ancestry_munge(x):
    checkpoint_output = checkpoints.ancestry_reporter.get(name=x, outdir=outdir).output[0]
    checkpoint_output = outdir + "/" + x + "/ancestry/ancestry_report.txt"
    ancestry_report_df = pd.read_table(checkpoint_output, sep=' ')
    return ancestry_report_df['population'].tolist()
    