#######
# Read target_list
#######

if 'target_list' in config:
  target_list_df = pd.read_table(config["target_list"], sep=r'\s+')
  target_list_df_23andMe = target_list_df.loc[target_list_df['type'] == '23andMe']
  samp_types = ['plink1', 'plink2', 'bgen', 'vcf']
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

# Check specific target paths exist
check_target_paths(df = target_list_df, chr = str(get_chr_range(config['testing'])[0]))

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
      lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'path'].iloc[0],
      rules.download_impute2_data.output
    output:
      "{outdir}/{name}/geno/imputed/{name}.chr{chr}.bed"
    conda:
      "../envs/analysis.yaml"
    params:
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
        --ref resources/data/impute2/1000GP_Phase3 \
        --n_core {resources.cpus}"

  rule impute_23andme_all_chr:
    input:
      lambda w: expand("{outdir}/{name}/geno/imputed/{name}.chr{chr}.bed", name=w.name, chr=get_chr_range(testing = config['testing']), outdir=outdir)
    output:
      touch("{outdir}/reference/target_checks/{name}/impute_23andme_all_chr.done")

  rule impute_23andme:
    input:
      lambda w: expand("{outdir}/reference/target_checks/{name}/impute_23andme_all_chr.done", name=target_list_df_23andMe['name'], outdir=outdir)

##
# Convert to PLINK and harmonise with reference
##

def format_target_input(outdir, name):
    inputs = []
    filtered_df = target_list_df.loc[target_list_df['name'] == name, 'type']
    if filtered_df.iloc[0] == '23andMe':
        inputs.append(f"{outdir}/reference/target_checks/{name}/impute_23andme_all_chr.done")
    return inputs

def target_prefix(outdir, name):
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == '23andMe':
      return outdir + "/" + name + "/geno/imputed/" + name
    else:
      path=target_list_df.loc[target_list_df['name'] == name, 'path'].iloc[0]
      return path

def target_path(outdir, name, chr):
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == '23andMe':
      return outdir + "/" + name + "/geno/imputed/" + name + ".chr" + chr + ".bed"
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == 'plink1':
      return target_list_df.loc[target_list_df['name'] == name, 'path'].iloc[0] + ".chr" + chr + ".bed"
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == 'plink2':
      return target_list_df.loc[target_list_df['name'] == name, 'path'].iloc[0] + ".chr" + chr + ".pgen"
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == 'bgen':
      return target_list_df.loc[target_list_df['name'] == name, 'path'].iloc[0] + ".chr" + chr + ".bgen"
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == 'vcf':
      return target_list_df.loc[target_list_df['name'] == name, 'path'].iloc[0] + ".chr" + chr + ".vcf.gz"

def target_type(name):
    if target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0] == '23andMe':
      return "plink1"
    else:
      return target_list_df.loc[target_list_df['name'] == name, 'type'].iloc[0]

if 'target_list' in config:
  rule format_target_i:
    input:
      lambda w: format_target_input(outdir = w.outdir, name = w.name),
      lambda w: target_path(outdir = w.outdir, name = w.name, chr = w.chr),
      rules.download_default_ref.output,
      rules.install_genoutils.output
    output:
      "{outdir}/{name}/geno/{name}.ref.chr{chr}.bed"
    conda:
      "../envs/analysis.yaml"
    params:
      config_file = config["config_file"],
      testing=config["testing"],
      prefix= lambda w: target_prefix(outdir = w.outdir, name = w.name),
      type= lambda w: target_type(name = w.name)
    shell:
      "Rscript ../Scripts/format_target/format_target.R \
        --target {params.prefix}.chr{wildcards.chr} \
        --format {params.type} \
        --ref resources/data/ref/ref.chr{wildcards.chr} \
        --output {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr{wildcards.chr}"

rule format_target_all_chr:
  input:
    lambda w: expand("{outdir}/{name}/geno/{name}.ref.chr{chr}.bed", name=w.name, chr=get_chr_range(testing = config['testing']), outdir=outdir)
  output:
    touch("{outdir}/reference/target_checks/{name}/format_target_all_chr.done")

rule format_target:
  input:
    lambda w: expand("{outdir}/reference/target_checks/{name}/format_target_all_chr.done", name=target_list_df['name'], outdir=outdir)

####
# Ancestry inference
####

rule ancestry_inference_i:
  input:
    "{outdir}/reference/target_checks/{name}/format_target_all_chr.done"
  output:
    touch("{outdir}/reference/target_checks/{name}/ancestry_inference.done")
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
  input: expand("{outdir}/reference/target_checks/{name}/ancestry_inference.done", name=target_list_df['name'], outdir=outdir)

# Create a file listing target samples and population assignments
checkpoint ancestry_reporter:
  input:
    "{outdir}/reference/target_checks/{name}/ancestry_inference.done",
  output:
    touch("{outdir}/reference/target_checks/{name}/ancestry_reporter.done")
  conda:
    "../envs/analysis.yaml"
  shell:
    "Rscript ../Scripts/pipeline_misc/ancestry_reporter.R {wildcards.name} {outdir}"

rule run_ancestry_reporter:
  input: expand("{outdir}/reference/target_checks/{name}/ancestry_reporter.done", name=target_list_df['name'], outdir=outdir)

####
# Population outlier detection and target sample specific PC calculation
####

rule outlier_detection_i:
  resources:
    mem_mb=15000
  input:
    "{outdir}/reference/target_checks/{name}/ancestry_reporter.done"
  output:
    touch('{outdir}/reference/target_checks/{name}/outlier_detection.done')
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
    lambda w: expand("{outdir}/reference/target_checks/{name}/outlier_detection.done", name=target_list_df_samp['name'], outdir=outdir)
