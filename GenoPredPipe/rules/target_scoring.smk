# Create a function summarising which populations are present in target      
def ancestry_munge(x):
    checkpoint_output = checkpoints.ancestry_reporter.get(name=x, outdir=outdir).output[0]
    checkpoint_output = outdir + "/" + x + "/ancestry/ancestry_report.txt"
    ancestry_report_df = pd.read_table(checkpoint_output, sep=' ')
    return ancestry_report_df['population'].tolist()
    
# Create a function that 
def get_score_file(w):
  if w.method != 'external':
    return f"{outdir}/resources/data/ref/pgs_score_files/{w.method}/{w.gwas}/ref.{w.gwas}.score.gz"
  else:
    # Assuming score_list_df is defined in the global scope
    return score_list_df.loc[score_list_df['name'] == w.gwas, 'path'].iloc[0]

####
# Projected PCs
####

rule pc_projection:
  input:
    "{outdir}/resources/data/target_checks/{name}/ancestry_reporter.done",
    rules.run_ref_pca.input,
    "../Scripts/target_scoring/target_scoring.R"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/pc_projection_{population}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/target_scoring/target_scoring.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --target_keep {outdir}/{wildcards.name}/ancestry/keep_files/model_based/{wildcards.population}.keep \
      --ref_freq_chr resources/data/ref/freq_files/{wildcards.population}/ref.{wildcards.population}.chr \
      --ref_score resources/data/ref/pc_score_files/{wildcards.population}/ref.{wildcards.population}.pcs.eigenvec.var.gz \
      --ref_scale resources/data/ref/pc_score_files/{wildcards.population}/ref.{wildcards.population}.pcs.{wildcards.population}.scale \
      --plink2 plink2 \
      --output {outdir}/{wildcards.name}/pcs/projected/{wildcards.population}/{wildcards.name}.{wildcards.population}"
      
rule run_pc_projection_all_pop:
  input: 
    lambda w: expand("{outdir}/resources/data/target_checks/{name}/pc_projection_{population}.done", name=w.name, population=ancestry_munge("{}".format(w.name)), outdir=outdir)
  output:
    touch("{outdir}/resources/data/target_checks/{name}/run_pc_projection_all_pop.done")

rule run_pc_projection_all:
  input: 
    expand("{outdir}/resources/data/target_checks/{name}/run_pc_projection_all_pop.done", name=target_list_df['name'], outdir=outdir)

####
# Polygenic scoring
####

rule target_pgs:
  input:
    "{outdir}/resources/data/target_checks/{name}/ancestry_reporter.done",
    "{outdir}/resources/data/ref/pgs_score_files/{method}/{gwas}/ref.{gwas}.EUR.scale",
    "../Scripts/target_scoring/target_scoring.R",
    "../Scripts/functions/misc.R"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/target_pgs_{method}_{population}_{gwas}.done")
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    score_file= lambda w: get_score_file(w)
  shell:
    "Rscript ../Scripts/target_scoring/target_scoring.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --target_keep {outdir}/{wildcards.name}/ancestry/keep_files/model_based/{wildcards.population}.keep \
      --ref_score {params.score_file} \
      --ref_scale {outdir}/resources/data/ref/pgs_score_files/{wildcards.method}/{wildcards.gwas}/ref.{wildcards.gwas}.{wildcards.population}.scale \
      --ref_freq_chr resources/data/ref/freq_files/{wildcards.population}/ref.{wildcards.population}.chr \
      --plink2 plink2 \
      --pheno_name {wildcards.gwas} \
      --output {outdir}/{wildcards.name}/pgs/{wildcards.population}/{wildcards.method}/{wildcards.gwas}/{wildcards.name}.{wildcards.gwas}.{wildcards.population}"

rule run_target_pgs_all_gwas:
  input:
    lambda w: expand("{outdir}/resources/data/target_checks/{name}/target_pgs_{method}_{population}_{gwas}.done", name=w.name, gwas= score_list_df['name'] if w.method == 'external' else (gwas_list_df['name'] if w.method in ['dbslmm','lassosum','megaprs'] else gwas_list_df['name']), population=w.population, method=w.method, outdir=outdir)
  output:
    touch("{outdir}/resources/data/target_checks/{name}/run_target_pgs_{method}_all_gwas_{population}.done")

rule run_target_pgs_all_pop:
  input: 
    lambda w: expand("{outdir}/resources/data/target_checks/{name}/run_target_pgs_{method}_all_gwas_{population}.done", name=w.name, population=ancestry_munge("{}".format(w.name)), method=w.method, outdir=outdir)
  output:
    touch("{outdir}/resources/data/target_checks/{name}/run_target_pgs_{method}_all_pop.done")

rule run_target_pgs_all_method:
  input: 
    lambda w: expand("{outdir}/resources/data/target_checks/{name}/run_target_pgs_{method}_all_pop.done", method=pgs_methods, name=w.name, outdir=outdir)
  output:
    touch("{outdir}/resources/data/target_checks/{name}/run_target_pgs_all_method.done")

rule run_target_pgs_all_name:
  input: 
    expand("{outdir}/resources/data/target_checks/{name}/run_target_pgs_all_method.done", name=target_list_df['name'], outdir=outdir)
  output:
    touch("{outdir}/resources/data/target_checks/run_target_pgs_all_name.done")
