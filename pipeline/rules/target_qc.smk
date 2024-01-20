#######
# Read target_list
#######

if 'target_list' in config:
  target_list_df = pd.read_table(config["target_list"], sep=r'\s+')
  target_list_df_23andMe = target_list_df.loc[target_list_df['type'] == '23andMe']
  samp_types = ['plink1', 'bgen', 'vcf']
  target_list_df_samp = target_list_df[target_list_df['type'].isin(samp_types)]
  target_list_df_indiv_report = target_list_df.loc[(target_list_df['indiv_report'].isin(['T', 'TRUE', True]))]
else:
  target_list_df = pd.DataFrame(columns = ["name", "path" "type", "indiv_report"])
  target_list_df_23andMe = pd.DataFrame(columns = ["name", "path" "type", "indiv_report"])
  target_list_df_samp = pd.DataFrame(columns = ["name", "path" "type", "indiv_report"])
  target_list_df_indiv_report = pd.DataFrame(columns = ["name", "path" "type", "indiv_report"])

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
      mem_mb=80000,
      cpus=config.get("ncores", 10),
      time_min=800
    input:
      config['target_list'],
      config['config_file'],
      rules.download_impute2_data.output,
      "../Scripts/23andMe_imputer/23andMe_imputer.R"
    output:
      touch("{outdir}/resources/data/target_checks/{name}/impute_23andme-{chr}.done")
    conda:
      "../envs/analysis.yaml"
    params:
      name= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'name'].iloc[0],
      path= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0]
    shell:
      "Rscript ../Scripts/23andMe_imputer/23andMe_imputer.R \
        --FID {params.name} \
        --IID {params.name} \
        --geno {params.path} \
        --output {outdir}/{params.name}/geno/imputed/{params.name}.chr{wildcards.chr} \
        --chr {wildcards.chr} \
        --ref resources/data/impute2/1000GP_Phase3 \
        --shapeit shapeit \
        --impute2 impute2 \
        --n_core {resources.cpus}"

  rule impute_23andme_all_chr:
    input:
      lambda w: expand("{outdir}/resources/data/target_checks/{name}/impute_23andme-{chr}.done", name=w.name, chr=get_chr_range(testing = config['testing']), outdir=outdir)
    output:
      touch("{outdir}/resources/data/target_checks/{name}/impute_23andme_all_chr.done")

  rule impute_23andme:
    input:
      lambda w: expand("{outdir}/resources/data/target_checks/{name}/impute_23andme_all_chr.done", name=target_list_df_23andMe['name'], outdir=outdir)

##
# Convert to PLINK and harmonise with reference
##

def format_target_input(outdir, name):
    inputs = []
    filtered_df = target_list_df.loc[target_list_df['name'] == name, 'type']
    if filtered_df.iloc[0] == '23andMe':
        inputs.append(f"{outdir}/resources/data/target_checks/{name}/impute_23andme_all_chr.done")
    return inputs

def target_path(outdir, name):
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == '23andMe':
      return outdir + "/" + name + "/geno/imputed/" + name
    else:
      path=target_list_df.loc[target_list_df['name'] == name, 'path'].iloc[0]
      return path

def target_type(name):
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == '23andMe':
      return "plink1"
    else:
      return target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0]

if 'target_list' in config:
  rule format_target_i:
    input:
      lambda w: format_target_input(outdir = w.outdir, name = w.name),
      config['target_list'],
      config['config_file'],
      rules.get_dependencies.output,
      "../Scripts/format_target/format_target.R"
    output:
      "{outdir}/{name}/geno/{name}.ref.chr{chr}.bed"
    conda:
      "../envs/analysis.yaml"
    params:
      path= lambda w: target_path(outdir = w.outdir, name = w.name),
      type= lambda w: target_type(name = w.name)
    shell:
      "Rscript ../Scripts/format_target/format_target.R \
        --target {params.path}.chr{wildcards.chr} \
        --format {params.type} \
        --ref resources/data/ref/ref.chr{wildcards.chr} \
        --output {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr{wildcards.chr}"

def get_chr_range(testing):
  if testing != 'NA':
    val = testing[-2:]
    val = int(val)
    return val
  else:
    return range(1, 23)  # Full range for normal operation

rule format_target_all_chr:
  input:
    lambda w: expand("{outdir}/{name}/geno/{name}.ref.chr{chr}.bed", name=w.name, chr=get_chr_range(testing = config['testing']), outdir=outdir)
  output:
    touch("{outdir}/resources/data/target_checks/{name}/format_target_all_chr.done")

rule format_target:
  input:
    lambda w: expand("{outdir}/resources/data/target_checks/{name}/format_target_all_name.done", name=target_list_df['name'], outdir=outdir)

####
# Ancestry inference
####

rule ancestry_inference_i:
  input:
    "{outdir}/resources/data/target_checks/{name}/format_target_all_chr.done",
    "../Scripts/Ancestry_identifier/Ancestry_identifier.R"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/ancestry_inference.done")
  conda:
    "../envs/analysis.yaml"
  params:
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/Ancestry_identifier/Ancestry_identifier.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --ref_plink_chr resources/data/ref/ref.chr \
      --output {outdir}/{wildcards.name}/ancestry/{wildcards.name}.Ancestry \
      --pop_data resources/data/ref/ref.pop.txt \
      --test {params.testing}"

rule ancestry_inference:
  input: expand("{outdir}/resources/data/target_checks/{name}/ancestry_inference.done", name=target_list_df['name'], outdir=outdir)

# Create a file listing target samples and population assignments
checkpoint ancestry_reporter:
  input:
    "{outdir}/resources/data/target_checks/{name}/ancestry_inference.done",
    "scripts/ancestry_reporter.R"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/ancestry_reporter.done")
  conda:
    "../envs/analysis.yaml"
  shell:
    "Rscript scripts/ancestry_reporter.R {wildcards.name} {outdir}"

rule run_ancestry_reporter:
  input: expand("{outdir}/resources/data/target_checks/{name}/ancestry_reporter.done", name=target_list_df['name'], outdir=outdir)

####
# Population outlier detection and target sample specific PC calculation
####

rule outlier_detection_i:
  resources:
    mem_mb=15000
  input:
    "{outdir}/resources/data/target_checks/{name}/ancestry_reporter.done",
    "../Scripts/outlier_detection/outlier_detection.R"
  output:
    touch('{outdir}/resources/data/target_checks/{name}/outlier_detection.done')
  conda:
    "../envs/analysis.yaml"
  params:
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/outlier_detection/outlier_detection.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --keep_list {outdir}/{wildcards.name}/ancestry/keep_list.txt \
      --test {params.testing} \
      --output {outdir}/{wildcards.name}/pcs/within_sample/{wildcards.name}.outlier_detection"

rule outlier_detection:
  input:
    lambda w: expand("{outdir}/resources/data/target_checks/{name}/outlier_detection.done", name=target_list_df_samp['name'], outdir=outdir)
