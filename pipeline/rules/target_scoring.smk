def ancestry_munge(x, scaling=['continuous'], target_populations=None):
    """
    Process ancestry report and return list of populations based on scaling type and target population filter.

    Parameters:
    - x (str): Name of the ancestry checkpoint
    - scaling (list): ['discrete', 'continuous'] specifying how to scale ancestry effects
    - target_populations (list or None): Optional list of population codes to include (e.g., ['EUR', 'EAS'])

    Returns:
    - list: A list of populations to be used
    """

    if not isinstance(scaling, list):
        raise ValueError("The scaling parameter must be a list (e.g., ['discrete', 'continuous']).")

    # Retrieve ancestry reporter output
    checkpoints.ancestry_reporter.get(name=x).output[0]
    checkpoint_output = outdir + "/" + x + "/ancestry/ancestry_report.txt"

    # Read ancestry report
    ancestry_report_df = pd.read_table(checkpoint_output, sep=' ')
    available_populations = ancestry_report_df['population'].tolist()
    
    # Extract population list
    population_list = ancestry_report_df['population'].tolist()

    # Filter based on target_populations if specified
    if target_populations not in [None, "NA"]:
        missing = [pop for pop in target_populations if pop not in available_populations]
        if missing:
            print(f"[ancestry_munge] The following target populations were not found in ancestry report: {missing}")
        population_list = [pop for pop in target_populations if pop in available_populations]
    else:
        population_list = available_populations

    # Handle scaling logic
    if 'continuous' in scaling and 'discrete' not in scaling:
        return ['TRANS']
    elif 'continuous' in scaling and 'discrete' in scaling:
        return population_list + ['TRANS']
    elif 'discrete' in scaling and 'continuous' not in scaling:
        return population_list
    else:
        raise ValueError("Invalid value for scaling. Must include 'continuous', 'discrete', or both.")

# Define which pgs_methods are can be applied to any GWAS population
pgs_methods_noneur = ['ptclump','lassosum','megaprs','prscs','dbslmm']

####
# Projected PCs
####

rule pc_projection_i:
  input:
    f"{outdir}/reference/target_checks/{{name}}/ancestry_reporter.done",
    f"{outdir}/reference/pc_score_files/{{population}}/ref-{{population}}-pcs.EUR.scale"
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
    target_keep=lambda wildcards: "NA" if wildcards.population == "TRANS" else f"{outdir}/{wildcards.name}/ancestry/keep_files/model_based/{wildcards.population}.keep",
    sbayesrc_ldref=lambda wildcards: f"{sbayesrc_ldref}/{wildcards.population}"
  shell:
    "Rscript ../Scripts/target_scoring/target_scoring.R \
      --target_plink_chr {outdir}/{wildcards.name}/geno/{wildcards.name}.ref.chr \
      --target_keep {params.target_keep} \
      --ref_freq_chr {refdir}/freq_files/{wildcards.population}/ref.{wildcards.population}.chr \
      --ref_score {outdir}/reference/pc_score_files/{wildcards.population}/ref-{wildcards.population}-pcs.eigenvec.var.gz \
      --ref_scale {outdir}/reference/pc_score_files/{wildcards.population}/ref-{wildcards.population}-pcs.{wildcards.population}.scale \
      --ref_ld_eigen {params.sbayesrc_ldref} \
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
    time_min=2800
  threads: config['cores_target_pgs']
  input:
    f"{outdir}/reference/target_checks/{{name}}/ancestry_reporter.done",
    lambda w: f"{outdir}/reference/target_checks/{{name}}/pc_projection-TRANS.done" if w.population == "TRANS" else [],
    rules.ref_pgs.output
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
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/target_pgs-{{population}}.done", name=w.name, population = ancestry_munge(w.name, scaling = config["pgs_scaling"], target_populations = config['target_populations']))
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/target_pgs.done")

rule target_pgs:
  input:
    expand(f"{outdir}/reference/target_checks/{{name}}/target_pgs.done", name=target_list_df['name'])

