def ancestry_munge(x, scaling='continuous'):
    # Ensure scaling is a list
    if not isinstance(scaling, list):
        raise ValueError("The scaling parameter must be a list (e.g., ['discrete', 'continuous']).")

    # Retrieve ancestry reporter output
    checkpoints.ancestry_reporter.get(name=x).output[0]
    checkpoint_output = outdir + "/" + x + "/ancestry/ancestry_report.txt"

    # Read ancestry report
    ancestry_report_df = pd.read_table(checkpoint_output, sep=' ')

    # Extract population list
    population_list = ancestry_report_df['population'].tolist()

    # Handle scaling logic
    if 'continuous' in scaling and 'discrete' not in scaling:
        # Only continuous scaling
        return ['TRANS']
    elif 'continuous' in scaling and 'discrete' in scaling:
        # Both continuous and discrete scaling
        population_list.append('TRANS')
        return population_list
    elif 'discrete' in scaling and 'continuous' not in scaling:
        # Only discrete scaling
        return population_list
    else:
        # Raise an error for invalid scaling input
        raise ValueError("Invalid value for scaling. Must include 'continuous', 'discrete', or both.")

# Define which pgs_methods are can be applied to any GWAS population
pgs_methods_noneur = ['ptclump','lassosum','megaprs','prscs','dbslmm']

####
# Projected PCs
####

rule pc_projection_i:
  input:
    f"{outdir}/reference/target_checks/{{name}}/ancestry_reporter.done",
    f"{resdir}/data/ref/pc_score_files/{{population}}/ref-{{population}}-pcs.EUR.scale"
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/pc_projection-{{population}}.done")
  benchmark:
    f"{outdir}/reference/benchmarks/pc_projection_i-{{name}}-{{population}}.txt"
  log:
    f"{outdir}/reference/logs/pc_projection_i-{{name}}-{{population}}.log"
  conda:
    "../envs/analysis.yaml"
  params:
    testing=config["testing"],
    target_keep=lambda wildcards: "NA" if wildcards.population == "TRANS" else f"{outdir}/{wildcards.name}/ancestry/keep_files/model_based/{wildcards.population}.keep"
  shell:
    "Rscript ../Scripts/target_scoring/target_scoring.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --target_keep {params.target_keep} \
      --ref_freq_chr {refdir}/freq_files/{wildcards.population}/ref.{wildcards.population}.chr \
      --ref_score {resdir}/data/ref/pc_score_files/{wildcards.population}/ref-{wildcards.population}-pcs.eigenvec.var.gz \
      --ref_scale {resdir}/data/ref/pc_score_files/{wildcards.population}/ref-{wildcards.population}-pcs.{wildcards.population}.scale \
      --plink2 plink2 \
      --test {params.testing} \
      --output {outdir}/{wildcards.name}/pcs/projected/{wildcards.population}/{wildcards.name}-{wildcards.population} > {log} 2>&1"

rule pc_projection_all:
  input:
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/pc_projection-{{population}}.done", name=w.name, population=ancestry_munge("{}".format(w.name), scaling = config["pgs_scaling"]))
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/pc_projection.done")

rule pc_projection:
  input:
    expand(f"{outdir}/reference/target_checks/{{name}}/pc_projection.done", name=target_list_df['name'])

####
# Polygenic scoring
####

rule target_pgs_i:
  resources:
    mem_mb=config['mem_target_pgs'],
    time_min=1600
  threads: config['cores_target_pgs']
  input:
    f"{outdir}/reference/target_checks/{{name}}/ancestry_reporter.done",
    lambda w: f"{outdir}/reference/target_checks/{{name}}/pc_projection-TRANS.done" if w.population == "TRANS" else [],
    rules.prep_pgs.input
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/target_pgs-{{population}}.done")
  benchmark:
    f"{outdir}/reference/benchmarks/target_pgs_i-{{name}}-{{population}}.txt"
  log:
    f"{outdir}/reference/logs/target_pgs_i-{{name}}-{{population}}.log"
  conda:
    "../envs/analysis.yaml"
  params:
    testing=config["testing"],
    config_file = config["config_file"]
  shell:
    "Rscript ../Scripts/target_scoring/target_scoring_pipeline.R \
      --config {params.config_file} \
      --name {wildcards.name} \
      --population {wildcards.population} \
      --plink2 plink2 \
      --test {params.testing} \
      --n_cores {threads} > {log} 2>&1"

rule target_pgs_all:
  input:
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/target_pgs-{{population}}.done", name=w.name, population = ancestry_munge(w.name, scaling = config["pgs_scaling"]))
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/target_pgs.done")

rule target_pgs:
  input:
    expand(f"{outdir}/reference/target_checks/{{name}}/target_pgs.done", name=target_list_df['name'])

