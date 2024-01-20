# Create a function summarising which populations are present in target
def ancestry_munge(x):
    checkpoint_output = checkpoints.ancestry_reporter.get(name=x, outdir=outdir).output[0]
    checkpoint_output = outdir + "/" + x + "/ancestry/ancestry_report.txt"
    ancestry_report_df = pd.read_table(checkpoint_output, sep=' ')
    return ancestry_report_df['population'].tolist()

# Create a function summarising which score files matched sufficiently with reference
def score_munge(outdir):
    checkpoint_output = checkpoints.score_reporter.get(outdir = outdir).output[0]
    checkpoint_output = outdir + "/resources/data/ref/pgs_score_files/external/score_report.txt"
    score_report_df = pd.read_table(checkpoint_output, sep=' ')
    score_report_df = score_report_df[(score_report_df['pass'].isin(['T', 'TRUE', True]))]
    return score_report_df['name'].tolist()

####
# Projected PCs
####

rule pc_projection_i:
  input:
    "{outdir}/resources/data/target_checks/{name}/ancestry_reporter.done",
    rules.ref_pca.input,
    "../Scripts/target_scoring/target_scoring.R"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/pc_projection-{population}.done")
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
      --output {outdir}/{wildcards.name}/pcs/projected/{wildcards.population}/{wildcards.name}-{wildcards.population}"

rule pc_projection_all_pop:
  input:
    lambda w: expand("{outdir}/resources/data/target_checks/{name}/pc_projection-{population}.done", name=w.name, population=ancestry_munge("{}".format(w.name)), outdir=outdir)
  output:
    touch("{outdir}/resources/data/target_checks/{name}/pc_projection_all_pop.done")

rule pc_projection:
  input:
    expand("{outdir}/resources/data/target_checks/{name}/pc_projection_all_pop.done", name=target_list_df['name'], outdir=outdir)

####
# Polygenic scoring
####

def check_pgs(w):
  # Check if the path value is not NA
  if w.method == 'external':
      return ["{outdir}/resources/data/target_checks/prep_pgs_external_i-{gwas}.done"]
  else:
      return ["{outdir}/resources/data/ref/pgs_score_files/{method}/{gwas}/ref-{gwas}-EUR.scale"]
      
rule target_pgs_i:
  input:
    "{outdir}/resources/data/target_checks/{name}/ancestry_reporter.done",
    lambda w: check_pgs(w),
    "../Scripts/target_scoring/target_scoring.R",
    "../Scripts/functions/misc.R"
  output:
    touch("{outdir}/resources/data/target_checks/{name}/target_pgs-{method}-{population}-{gwas}.done")
  conda:
    "../envs/analysis.yaml"
  params:
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/target_scoring/target_scoring.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --target_keep {outdir}/{wildcards.name}/ancestry/keep_files/model_based/{wildcards.population}.keep \
      --ref_score {outdir}/resources/data/ref/pgs_score_files/{wildcards.method}/{wildcards.gwas}/ref-{wildcards.gwas}.score.gz \
      --ref_scale {outdir}/resources/data/ref/pgs_score_files/{wildcards.method}/{wildcards.gwas}/ref-{wildcards.gwas}-{wildcards.population}.scale \
      --ref_freq_chr resources/data/ref/freq_files/{wildcards.population}/ref.{wildcards.population}.chr \
      --plink2 plink2 \
      --pheno_name {wildcards.gwas} \
      --test {params.testing} \
      --output {outdir}/{wildcards.name}/pgs/{wildcards.population}/{wildcards.method}/{wildcards.gwas}/{wildcards.name}-{wildcards.gwas}-{wildcards.population}"

rule target_pgs_all_gwas:
  input:
    lambda w: expand("{outdir}/resources/data/target_checks/{name}/target_pgs-{method}-{population}-{gwas}.done", name=w.name, gwas=score_munge(outdir = outdir) if w.method == 'external' else (gwas_list_df['name'] if w.method in ['ptclump','lassosum','megaprs'] else gwas_list_df_eur['name']), population=w.population, method=w.method, outdir=outdir)
  output:
    touch("{outdir}/resources/data/target_checks/{name}/target_pgs-{method}-{population}.done")

rule target_pgs_all_pop:
  input:
    lambda w: expand("{outdir}/resources/data/target_checks/{name}/target_pgs-{method}-{population}.done", name=w.name, population=ancestry_munge("{}".format(w.name)), method=w.method, outdir=outdir)
  output:
    touch("{outdir}/resources/data/target_checks/{name}/target_pgs-{method}.done")

rule target_pgs_all_method:
  input:
    lambda w: expand("{outdir}/resources/data/target_checks/{name}/target_pgs-{method}.done", method=pgs_methods, name=w.name, outdir=outdir)
  output:
    touch("{outdir}/resources/data/target_checks/{name}/target_pgs.done")

rule target_pgs:
  input:
    expand("{outdir}/resources/data/target_checks/{name}/target_pgs.done", name=target_list_df['name'], outdir=outdir)