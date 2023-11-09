####
# PC scoring
####

rule target_pc:
  input:
    "resources/data/target_checks/{name}/ancestry_reporter.done",
    rules.run_pop_pc_scoring.input,
    "../Scripts/target_scoring/target_scoring.R"
  output:
    touch("resources/data/target_checks/{name}/target_pc_{population}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0]
  shell:
    "Rscript ../Scripts/target_scoring/target_scoring.R \
      --target_plink_chr {params.output}/{wildcards.name}/{wildcards.name}.ref.chr \
      --target_keep {params.output}/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.{wildcards.population}.keep \
      --ref_freq_chr resources/data/ref/freq_files/{wildcards.population}/ref.{wildcards.population}.chr \
      --ref_score resources/data/ref/pc_score_files/{wildcards.population}/ref.{wildcards.population}.eigenvec.var \
      --ref_scale resources/data/ref/pc_score_files/{wildcards.population}/ref.{wildcards.population}.scale \
      --plink2 plink2 \
      --output {params.output}/{wildcards.name}/projected_pc/{wildcards.population}/{wildcards.name}"
      
rule run_target_pc_all_pop:
  input: 
    lambda w: expand("resources/data/target_checks/{name}/target_pc_{population}.done", name=w.name, population=ancestry_munge("{}".format(w.name)))
  output:
    touch("resources/data/target_checks/{name}/run_target_pc_all_pop.done")

rule run_target_pc_all:
  input: 
    expand("resources/data/target_checks/{name}/run_target_pc_all_pop.done", name=target_list_df['name'])

####
# Polygenic scoring
####

if score_list_file.is_file():
  pgs_methods = config['pgs_methods']  # Retrieve initial methods from config
  pgs_methods.append('external')
else:
  pgs_methods = config['pgs_methods']  # Retrieve initial methods from config

def get_score_file(w):
  if w.method != 'external':
    return f"resources/data/ref/prs_score_files/{w.method}/{w.gwas}/ref.{w.gwas}.score.gz"
  else:
    # Assuming score_list_df is defined in the global scope
    return score_list_df.loc[score_list_df['name'] == w.gwas, 'path'].iloc[0]


rule target_prs:
  input:
    "resources/data/target_checks/{name}/ancestry_reporter.done",
    "resources/data/ref/prs_score_files/{method}/{gwas}/ref.{gwas}.EUR.scale",
    "../Scripts/target_scoring/target_scoring.R",
    "../Scripts/functions/misc.R"
  output:
    touch("resources/data/target_checks/{name}/target_prs_{method}_{population}_{gwas}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    score_file= lambda w: get_score_file(w),
    output= lambda w: target_list_df.loc[target_list_df['name'] == "{}".format(w.name), 'output'].iloc[0]
  shell:
    "Rscript ../Scripts/target_scoring/target_scoring.R \
      --target_plink_chr {params.output}/{wildcards.name}/{wildcards.name}.ref.chr \
      --target_keep {params.output}/{wildcards.name}/ancestry/ancestry_all/{wildcards.name}.Ancestry.model_pred.{wildcards.population}.keep \
      --ref_score {params.score_file} \
      --ref_scale resources/data/ref/prs_score_files/{wildcards.method}/{wildcards.gwas}/ref.{wildcards.gwas}.{wildcards.population}.scale \
      --ref_freq_chr resources/data/ref/freq_files/{wildcards.population}/ref.{wildcards.population}.chr \
      --plink2 plink2 \
      --pheno_name {wildcards.gwas} \
      --output {params.output}/{wildcards.name}/prs/{wildcards.population}/{wildcards.method}/{wildcards.gwas}/{wildcards.name}.{wildcards.gwas}.{wildcards.population}"

rule run_target_prs_all_gwas:
  input:
    config["gwas_list"],
    lambda w: expand("resources/data/target_checks/{name}/target_prs_{method}_{population}_{gwas}.done", name=w.name, gwas= score_list_df['name'] if w.method == 'external' else (gwas_list_df['name'] if w.method in ['dbslmm','lassosum','megaprs'] else gwas_list_df['name']), population=w.population, method=w.method)
  output:
    touch("resources/data/target_checks/{name}/run_target_prs_{method}_all_gwas_{population}.done")

rule run_target_prs_all_pop:
  input: 
    lambda w: expand("resources/data/target_checks/{name}/run_target_prs_{method}_all_gwas_{population}.done", name=w.name, population=ancestry_munge("{}".format(w.name)), method=w.method)
  output:
    touch("resources/data/target_checks/{name}/run_target_prs_{method}_all_pop.done")

rule run_target_prs_all_method:
  input: 
    lambda w: expand("resources/data/target_checks/{name}/run_target_prs_{method}_all_pop.done", method=pgs_methods, name=w.name)
  output:
    touch("resources/data/target_checks/{name}/run_target_prs_all_method.done")

rule run_target_prs_all_name:
  input: 
    expand("resources/data/target_checks/{name}/run_target_prs_all_method.done", name=target_list_df['name'])
  output:
    touch("resources/data/target_checks/run_target_prs_all_name.done")
