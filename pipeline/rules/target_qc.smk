#######
# Read target_list
#######

####
# Format target data
####

##
# 23andMe
##
# Largely based on Impute.Me by Lasse Folkersen

if 'target_list' in config:
  rule impute_23andme_i:
    resources:
      mem_mb=10000*config['cores_impute_23andme'],
      time_min=800
    threads: config['cores_impute_23andme']
    input:
      lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0],
      rules.download_impute2_data.output,
      f"{resdir}/last_version.txt"
    output:
      f"{outdir}/{{name}}/geno/imputed/{{name}}.chr{{chr}}.bed"
    benchmark:
      f"{outdir}/reference/benchmarks/impute_23andme_i-{{name}}-{{chr}}.txt"
    log:
      f"{outdir}/reference/logs/impute_23andme_i-{{name}}-{{chr}}.log"
    conda:
      "../envs/analysis.yaml"
    params:
      outdir=config["outdir"],
      config_file = config["config_file"],
      name= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'name'].iloc[0],
      path= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0]
    shell:
      "Rscript ../Scripts/23andMe_imputer/23andMe_imputer.R \
        --FID {params.name} \
        --IID {params.name} \
        --geno {params.path} \
        --output {outdir}/{params.name}/geno/imputed/{params.name}.chr{wildcards.chr} \
        --chr {wildcards.chr} \
        --ref {resdir}/data/impute2/1000GP_Phase3 \
        --n_core {threads} > {log} 2>&1"

  rule impute_23andme_all:
    input:
      lambda w: expand(f"{outdir}/{{name}}/geno/imputed/{{name}}.chr{{chr}}.bed", name=w.name, chr=get_chr_range(testing = config['testing']))
    output:
      touch(f"{outdir}/reference/target_checks/{{name}}/impute_23andme.done")

  rule impute_23andme:
    input:
      lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/impute_23andme.done", name=target_list_df_23andMe['name'])

##
# Convert to PLINK and harmonise with reference
##

def format_target_input(name):
    inputs = []
    filtered_df = target_list_df.loc[target_list_df['name'] == name, 'type']
    if filtered_df.iloc[0] == '23andMe':
        inputs.append(f"{outdir}/reference/target_checks/{name}/impute_23andme.done")
    return inputs

def target_prefix(name):
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == '23andMe':
      return f"{outdir}/{name}/geno/imputed/{name}"
    else:
      path=target_list_df.loc[target_list_df['name'] == name, 'path'].iloc[0]
      return path

def target_path(name, chr):
    path=target_list_df.loc[target_list_df['name'] == name, 'path'].iloc[0]
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == '23andMe':
      return f"{outdir}/{name}/geno/imputed/{name}.chr{chr}.bed"
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == 'plink1':
      return f"{path}.chr{chr}.bed"
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == 'plink2':
      return f"{path}.chr{chr}.pgen"
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == 'bgen':
      return f"{path}.chr{chr}.bgen"
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == 'vcf':
      return f"{path}.chr{chr}.vcf.gz"

def target_type(name):
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == '23andMe':
      return "plink1"
    else:
      return target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0]

if 'target_list' in config:
  rule format_target_i:
    input:
      lambda w: format_target_input(name = w.name),
      lambda w: target_path(name = w.name, chr = w.chr),
      ref_input,
      rules.install_genoutils.output,
      f"{resdir}/last_version.txt"
    output:
      f"{outdir}/{{name}}/geno/{{name}}.ref.chr{{chr}}.pgen"
    benchmark:
      f"{outdir}/reference/benchmarks/format_target_i-{{name}}-{{chr}}.txt"
    log:
      f"{outdir}/reference/logs/format_target_i-{{name}}-{{chr}}.log"
    conda:
      "../envs/analysis.yaml"
    params:
      outdir=config["outdir"],
      refdir=config["refdir"],
      config_file = config["config_file"],
      testing=config["testing"],
      prefix= lambda w: target_prefix(name = w.name),
      type= lambda w: target_type(name = w.name)
    shell:
      "Rscript ../Scripts/format_target/format_target.R \
        --target {params.prefix}.chr{wildcards.chr} \
        --format {params.type} \
        --ref {refdir}/ref.chr{wildcards.chr} \
        --output {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr{wildcards.chr} > {log} 2>&1"

rule format_target_all:
  input:
    lambda w: expand(f"{outdir}/{{name}}/geno/{{name}}.ref.chr{{chr}}.pgen", name=w.name, chr=get_chr_range(testing = config['testing']))
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/format_target.done")

rule format_target:
  input:
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/format_target.done", name=target_list_df['name'])

####
# Ancestry inference
####

rule ancestry_inference_i:
  input:
    f"{outdir}/reference/target_checks/{{name}}/format_target.done"
  output:
    f"{outdir}/{{name}}/ancestry/{{name}}.Ancestry.model_pred"
  benchmark:
    f"{outdir}/reference/benchmarks/ancestry_inference_i-{{name}}.txt"
  log:
    f"{outdir}/reference/logs/ancestry_inference_i-{{name}}.log"
  conda:
    "../envs/analysis.yaml"
  params:
    testing=config["testing"],
    ancestry_prob_thresh=config["ancestry_prob_thresh"]
  shell:
    "rm -r -f {outdir}/{wildcards.name}/ancestry; \
     rm -f {outdir}/reference/target_checks/{wildcards.name}/ancestry_reporter.done; \
     Rscript ../Scripts/Ancestry_identifier/Ancestry_identifier.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --ref_plink_chr {refdir}/ref.chr \
      --output {outdir}/{wildcards.name}/ancestry/{wildcards.name}.Ancestry \
      --pop_data {refdir}/ref.pop.txt \
      --prob_thresh {params.ancestry_prob_thresh} \
      --test {params.testing} > {log} 2>&1"

rule ancestry_inference:
  input: expand(f"{outdir}/{{name}}/ancestry/{{name}}.Ancestry.model_pred", name=target_list_df['name'])

# Create a file listing target samples and population assignments
checkpoint ancestry_reporter:
  input:
    f"{outdir}/{{name}}/ancestry/{{name}}.Ancestry.model_pred",
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/ancestry_reporter.done")
  benchmark:
    f"{outdir}/reference/benchmarks/ancestry_reporter-{{name}}.txt"
  log:
    f"{outdir}/reference/logs/ancestry_reporter-{{name}}.log"
  conda:
    "../envs/analysis.yaml"
  shell:
    "Rscript ../Scripts/pipeline_misc/ancestry_reporter.R {wildcards.name} {outdir} > {log} 2>&1"

rule run_ancestry_reporter:
  input: expand(f"{outdir}/reference/target_checks/{{name}}/ancestry_reporter.done", name=target_list_df['name'])

####
# Population outlier detection and target sample specific PC calculation
####

rule outlier_detection_i:
  threads: config['cores_outlier_detection']
  resources:
    mem_mb=15000
  input:
    f"{outdir}/reference/target_checks/{{name}}/ancestry_reporter.done"
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/outlier_detection.done")
  benchmark:
    f"{outdir}/reference/benchmarks/outlier_detection_i-{{name}}.txt"
  log:
    f"{outdir}/reference/logs/outlier_detection_i-{{name}}.log"
  conda:
    "../envs/analysis.yaml"
  params:
    testing=config["testing"],
    unrel= lambda w: target_list_df.loc[target_list_df['name'] == w.name, 'unrel'].iloc[0]
  shell:
    "Rscript ../Scripts/outlier_detection/outlier_detection.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --keep_list {outdir}/{wildcards.name}/ancestry/keep_list.txt \
      --unrel {params.unrel} \
      --n_cores {threads} \
      --test {params.testing} \
      --output {outdir}/{wildcards.name}/pcs/within_sample/{wildcards.name}.outlier_detection > {log} 2>&1"

rule outlier_detection:
  input:
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/outlier_detection.done", name=target_list_df_samp['name'])
