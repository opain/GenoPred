import pandas as pd
target_list_df = pd.read_table(config["target_list"], sep=r'\s+')

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
# Formatting without imputation
##
# Here we will use a script that can accept multiple genetic data formats, insert RSIDs if not present, extract HapMap3 SNPs and output plink hard-call binaries.

rule format_target:
  input:
    config['target_list'],
    rules.get_dependencies.output,
    "../Scripts/geno_to_plink/geno_to_plink.R"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/format_target_{chr}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    path= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0],
    type= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'type'].iloc[0]
  shell:
    "Rscript ../Scripts/geno_to_plink/geno_to_plink.R \
      --target {params.path}.chr{wildcards.chr} \
      --format {params.type} \
      --ref resources/data/ref/ref.chr{wildcards.chr} \
      --plink2 plink2 \
      --output {outdir}/{wildcards.name}/geno/temp/{wildcards.name}.hm3.chr{wildcards.chr}"

rule run_format_target:
  input: 
    lambda w: expand("{outdir}/resources/data/target_checks/{name}/format_target_{chr}.done", name=w.name, chr=range(1, 23), outdir=outdir)
  output:
    touch("{outdir}/resources/data/target_checks/{name}/format_target.done")

rule run_format_target_2:
  input: 
    expand("{outdir}/resources/data/target_checks/{name}/format_target.done", name=target_list_df_samp_imp['name'], outdir=outdir)

####
# Harmonise with reference
####

rule harmonise_target_with_ref:
  input:
    rules.get_dependencies.output,
    lambda w: outdir + "/resources/data/target_checks/" + w.name + "/format_impute_23andme_target.done" if target_list_df.loc[target_list_df['name'] == w.name, 'type'].iloc[0] == '23andMe' else (outdir + "/resources/data/target_checks/{name}/format_target.done" if target_list_df.loc[target_list_df['name'] == w.name, 'type'].iloc[0] == 'samp_imp_plink1' else (outdir + "/resources/data/target_checks/{name}/format_target.done" if target_list_df.loc[target_list_df['name'] == w.name, 'type'].iloc[0] == 'samp_imp_bgen' else outdir + "/resources/data/target_checks/{name}/format_target.done")),
    "../Scripts/hm3_harmoniser/hm3_harmoniser.R"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/harmonise_target_with_ref.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/hm3_harmoniser/hm3_harmoniser.R \
      --target {outdir}/{wildcards.name}/geno/{wildcards.name}.hm3.chr \
      --ref resources/data/ref/ref.chr \
      --plink plink \
      --out {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr"

rule run_harmonise_target_with_ref:
  input: expand("{outdir}/resources/data/target_checks/{name}/harmonise_target_with_ref.done", name=target_list_df['name'], outdir=outdir)
  
# Delete temporary files
rule delete_temp_target_samp_files:
  input:
    "{outdir}/resources/data/target_checks/{name}/harmonise_target_with_ref.done"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/delete_temp_target_samp_files.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "rm {outdir}/{wildcards.name}/geno/{wildcards.name}.hm3.chr*"

####
# Population classification
####

rule target_population:
  input:
    "{outdir}/resources/data/target_checks/{name}/harmonise_target_with_ref.done",
    lambda w: outdir + "/resources/data/target_checks/" + w.name + "/harmonise_target_with_ref.done" if target_list_df.loc[target_list_df['name'] == w.name, 'type'].iloc[0] == '23andMe' else outdir + "/resources/data/target_checks/" + w.name + "/delete_temp_target_samp_files.done",
    "../Scripts/Ancestry_identifier/Ancestry_identifier.R"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/target_population.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/Ancestry_identifier/Ancestry_identifier.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --ref_plink_chr resources/data/ref/ref.chr \
      --n_pcs 6 \
      --plink plink \
      --plink2 plink2 \
      --output {outdir}/{wildcards.name}/ancestry/{wildcards.name}.Ancestry \
      --ref_pop_scale resources/data/ref/ref.keep.list \
      --pop_data resources/data/ref/ref.pop.txt \
      --prob_thresh 0.95"

rule run_target_population:
  input: expand("{outdir}/resources/data/target_checks/{name}/target_population.done", name=target_list_df['name'], outdir=outdir)

# Create a file listing target samples and population assignments
checkpoint ancestry_reporter:
  input:
    "{outdir}/resources/data/target_checks/{name}/target_population.done",
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

rule target_population_outlier_detection:
  resources: 
    mem_mb=15000
  input:
    "{outdir}/resources/data/target_checks/{name}/ancestry_reporter.done",
    "../Scripts/Population_outlier/Population_outlier.R"
  output:
    touch('{outdir}/resources/data/target_checks/{name}/target_population_outlier_detection.done')
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "ls {outdir}/{wildcards.name}/ancestry/{wildcards.name}.Ancestry.model_pred.*.keep > {outdir}/{wildcards.name}/ancestry/{wildcards.name}.Ancestry.model_pred.keep_list; Rscript ../Scripts/Population_outlier/Population_outlier.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --target_keep {outdir}/{wildcards.name}/ancestry/{wildcards.name}.Ancestry.model_pred.keep_list \
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
    