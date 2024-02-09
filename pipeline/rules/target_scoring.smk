# Create a function summarising which populations are present in target
def ancestry_munge(x):
    checkpoint_output = checkpoints.ancestry_reporter.get(name=x).output[0]
    checkpoint_output = outdir + "/" + x + "/ancestry/ancestry_report.txt"
    ancestry_report_df = pd.read_table(checkpoint_output, sep=' ')
    return ancestry_report_df['population'].tolist()

# Create a function summarising which score files matched sufficiently with reference
def score_munge():
    checkpoint_output = checkpoints.score_reporter.get().output[0]
    checkpoint_output = outdir + "/reference/pgs_score_files/external/score_report.txt"

    if os.path.exists(checkpoint_output):
      score_report_df = pd.read_table(checkpoint_output, sep=' ')
      score_report_df = score_report_df[(score_report_df['pass'].isin(['T', 'TRUE', True]))]
    else:
      score_report_df = pd.DataFrame(columns=['name', 'pass'])

    return score_report_df['name'].tolist()

####
# Projected PCs
####

rule pc_projection_i:
  input:
    f"{outdir}/reference/target_checks/{{name}}/ancestry_reporter.done",
    "resources/data/ref/pc_score_files/{population}/ref-{population}-pcs.EUR.scale"
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/pc_projection-{{population}}.done")
  benchmark:
    f"{outdir}/reference/benchmarks/pc_projection_i-{{name}}-{{population}}.txt"
  log:
    f"{outdir}/reference/logs/pc_projection_i-{{name}}-{{population}}.log"
  conda:
    "../envs/analysis.yaml"
  shell:
    "Rscript ../Scripts/target_scoring/target_scoring.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --target_keep {outdir}/{wildcards.name}/ancestry/keep_files/model_based/{wildcards.population}.keep \
      --ref_freq_chr resources/data/ref/freq_files/{wildcards.population}/ref.{wildcards.population}.chr \
      --ref_score resources/data/ref/pc_score_files/{wildcards.population}/ref-{wildcards.population}-pcs.eigenvec.var.gz \
      --ref_scale resources/data/ref/pc_score_files/{wildcards.population}/ref-{wildcards.population}-pcs.{wildcards.population}.scale \
      --plink2 plink2 \
      --output {outdir}/{wildcards.name}/pcs/projected/{wildcards.population}/{wildcards.name}-{wildcards.population} > {log} 2>&1"

rule pc_projection_all_pop:
  input:
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/pc_projection-{{population}}.done", name=w.name, population=ancestry_munge("{}".format(w.name)))
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/pc_projection_all_pop.done")

rule pc_projection:
  input:
    expand(f"{outdir}/reference/target_checks/{{name}}/pc_projection_all_pop.done", name=target_list_df['name'])

####
# Polygenic scoring
####

rule target_pgs_i:
  resources:
    mem_mb=8000,
    time_min=800
  threads: lambda w: 1 if w.method in ['dbslmm', 'sbayesr'] else 5
  input:
    f"{outdir}/reference/target_checks/{{name}}/ancestry_reporter.done",
    f"{outdir}/reference/target_checks/score_reporter.done",
    f"{outdir}/reference/target_checks/prep_pgs_{{method}}_i-{{gwas}}.done"
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/target_pgs-{{method}}-{{population}}-{{gwas}}.done")
  benchmark:
    f"{outdir}/reference/benchmarks/target_pgs_i-{{name}}-{{method}}-{{population}}-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/target_pgs_i-{{name}}-{{method}}-{{population}}-{{gwas}}.log"
  conda:
    "../envs/analysis.yaml"
  params:
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/target_scoring/target_scoring.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --target_keep {outdir}/{wildcards.name}/ancestry/keep_files/model_based/{wildcards.population}.keep \
      --ref_score {outdir}/reference/pgs_score_files/{wildcards.method}/{wildcards.gwas}/ref-{wildcards.gwas}.score.gz \
      --ref_scale {outdir}/reference/pgs_score_files/{wildcards.method}/{wildcards.gwas}/ref-{wildcards.gwas}-{wildcards.population}.scale \
      --ref_freq_chr resources/data/ref/freq_files/{wildcards.population}/ref.{wildcards.population}.chr \
      --plink2 plink2 \
      --pheno_name {wildcards.gwas} \
      --test {params.testing} \
      --n_cores {threads} \
      --output {outdir}/{wildcards.name}/pgs/{wildcards.population}/{wildcards.method}/{wildcards.gwas}/{wildcards.name}-{wildcards.gwas}-{wildcards.population} > {log} 2>&1"

rule target_pgs_all_gwas:
  input:
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/target_pgs-{{method}}-{{population}}-{{gwas}}.done", name=w.name, gwas = gwas_list_df['name'] if w.method in ['ptclump','lassosum','megaprs','prscs'] else gwas_list_df_eur['name'], population=w.population, method=w.method)
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/target_pgs_all_gwas-{{method}}-{{population}}.done")

rule target_pgs_all_method:
  input:
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/target_pgs_all_gwas-{{method}}-{{population}}.done", name=w.name, population=ancestry_munge("{}".format(w.name)), method=pgs_methods)
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/target_pgs_all_gwas-{{population}}.done")

rule target_pgs_all_score:
  input:
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/target_pgs-external-{{population}}-{{score}}.done", name=w.name, score=score_munge(), population=w.population)
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/target_pgs_all_score-{{population}}.done")

rule target_pgs_all_pop:
  input:
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/target_pgs_all_gwas-{{population}}.done", name=w.name, population=ancestry_munge("{}".format(w.name))),
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/target_pgs_all_score-{{population}}.done", name=w.name, population=ancestry_munge("{}".format(w.name)))
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/target_pgs.done")

rule target_pgs:
  input:
    expand(f"{outdir}/reference/target_checks/{{name}}/target_pgs.done", name=target_list_df['name'])

