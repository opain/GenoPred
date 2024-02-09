# Create PC score files specific to each population
rule ref_pca_i:
  input:
    rules.download_default_ref.output,
    rules.install_genoutils.output
  output:
    "resources/data/ref/pc_score_files/{population}/ref-{population}-pcs.EUR.scale"
  conda:
    "../envs/analysis.yaml",
  params:
    testing=config["testing"],
    outdir=config["outdir"]
  benchmark:
    "resources/data/benchmarks/ref_pca_i-{population}.txt"
  log:
    "resources/data/logs/ref_pca_i-{population}.log"
  shell:
    "Rscript ../Scripts/ref_pca/ref_pca.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{wildcards.population}.keep \
      --pop_data resources/data/ref/ref.pop.txt \
      --output resources/data/ref/pc_score_files/{wildcards.population}/ref-{wildcards.population}-pcs \
      --test {params.testing} > {log} 2>&1"

populations=["AFR","AMR","EAS","EUR","SAS"]

rule ref_pca:
  input: expand("resources/data/ref/pc_score_files/{population}/ref-{population}-pcs.EUR.scale", population=populations)

##
# QC and format GWAS summary statistics
##

# Read in the gwas_list or make an empty version
if 'gwas_list' in config and config["gwas_list"] != 'NA':
  gwas_list_df = pd.read_table(config["gwas_list"], sep=r'\s+')
else:
  gwas_list_df = pd.DataFrame(columns = ["name", "path", "population", "n", "sampling", "prevalence", "mean", "sd", "label"])

# Check whether gwas_list paths exist
check_list_paths(gwas_list_df)

# Identify gwas_list with population == 'EUR'
gwas_list_df_eur = gwas_list_df.loc[gwas_list_df['population'] == 'EUR']

if 'gwas_list' in config:
  rule sumstat_prep_i:
    input:
      rules.download_default_ref.output,
      rules.install_genoutils.output,
      lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0]
    output:
      f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz"
    benchmark:
      f"{outdir}/reference/benchmarks/sumstat_prep_i-{{gwas}}.txt"
    log:
      f"{outdir}/reference/logs/sumstat_prep_i-{{gwas}}.log"
    conda:
      "../envs/analysis.yaml"
    params:
      outdir=config["outdir"],
      config_file = config["config_file"],
      testing = config['testing'],
      population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
      n= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'n'].iloc[0],
      path= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0]
    shell:
      """
      sumstat_cleaner_script=$(Rscript -e 'cat(system.file("scripts", "sumstat_cleaner.R", package = "GenoUtils"))')
      Rscript $sumstat_cleaner_script \
        --sumstats {params.path} \
        --n {params.n} \
        --ref_chr resources/data/ref/ref.chr \
        --population {params.population} \
        --output {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned \
        --test {params.testing} > {log} 2>&1
      """

rule sumstat_prep:
  input: expand(f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz", gwas=gwas_list_df['name'])

##
# pT+clump (sparse, nested)
##

rule prep_pgs_ptclump_i:
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz"
  output:
    touch(f"{outdir}/reference/target_checks/prep_pgs_ptclump_i-{{gwas}}.done")
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_ptclump_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_ptclump_i-{{gwas}}.log"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/ptclump.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{params.population}.keep \
      --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --output {outdir}/reference/pgs_score_files/ptclump/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data resources/data/ref/ref.pop.txt \
      --test {params.testing} > {log} 2>&1"

rule prep_pgs_ptclump:
  input: expand(f"{outdir}/reference/target_checks/prep_pgs_ptclump_i-{{gwas}}.done", gwas=gwas_list_df['name'])

##
# DBSLMM
##

rule prep_pgs_dbslmm_i:
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.download_plink.output,
    rules.download_ldsc.output,
    rules.download_ldsc_ref.output,
    rules.download_hm3_snplist.output,
    rules.download_dbslmm.output,
    rules.download_ld_blocks.output
  output:
    touch(f"{outdir}/reference/target_checks/prep_pgs_dbslmm_i-{{gwas}}.done")
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_dbslmm_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_dbslmm_i-{{gwas}}.log"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    sampling= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'sampling'].iloc[0],
    prevalence= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'prevalence'].iloc[0],
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/dbslmm.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{params.population}.keep \
      --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --ld_blocks resources/data/ld_blocks/{params.population} \
      --plink resources/software/plink/plink \
      --dbslmm resources/software/dbslmm/software \
      --munge_sumstats resources/software/ldsc/munge_sumstats.py \
      --ldsc resources/software/ldsc/ldsc.py \
      --ldsc_ref resources/data/ldsc_ref/eur_w_ld_chr \
      --hm3_snplist resources/data/hm3_snplist/w_hm3.snplist \
      --sample_prev {params.sampling} \
      --pop_prev {params.prevalence} \
      --output {outdir}/reference/pgs_score_files/dbslmm/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data resources/data/ref/ref.pop.txt \
      --test {params.testing} > {log} 2>&1"

rule prep_pgs_dbslmm:
  input: expand(f"{outdir}/reference/target_checks/prep_pgs_dbslmm_i-{{gwas}}.done", gwas=gwas_list_df_eur['name'])

##
# PRScs
##
# Note. Threads are set to 1, and phi and chr are run in parallel. Increasing number of threads shows no improvement in speed.

# Set default values
n_cores_prscs = config.get("ncores", 10)

# Modify if the 'testing' condition is met
if config["testing"] != 'NA':
    n_cores_prscs = config.get("ncores", 5)

mem_prscs = 2000*n_cores_prscs

rule prep_pgs_prscs_i:
  resources:
    mem_mb=mem_prscs,
    time_min=800
  threads: n_cores_prscs
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.download_prscs_software.output,
    lambda w: "resources/data/prscs_ref/ldblk_ukbb_" + gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0].lower() + "/ldblk_ukbb_chr1.hdf5"
  output:
    touch(f"{outdir}/reference/target_checks/prep_pgs_prscs_i-{{gwas}}.done")
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_prscs_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_prscs_i-{{gwas}}.log"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0].lower(),
    testing=config["testing"]
  shell:
    """
    export MKL_NUM_THREADS=1; \
    export NUMEXPR_NUM_THREADS=1; \
    export OMP_NUM_THREADS=1; \
    export OPENBLAS_NUM_THREADS=1; \
    Rscript ../Scripts/pgs_methods/prscs.R \
    --ref_plink_chr resources/data/ref/ref.chr \
    --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
    --output {outdir}/reference/pgs_score_files/prscs/{wildcards.gwas}/ref-{wildcards.gwas} \
    --pop_data resources/data/ref/ref.pop.txt \
    --PRScs_path resources/software/prscs/PRScs.py \
    --PRScs_ref_path resources/data/prscs_ref/ldblk_ukbb_{params.population} \
    --n_cores {threads} \
    --phi_param 1e-6,1e-4,1e-2,1,auto \
    --test {params.testing} > {log} 2>&1
    """

rule prep_pgs_prscs:
  input: expand(f"{outdir}/reference/target_checks/prep_pgs_prscs_i-{{gwas}}.done", gwas=gwas_list_df['name'])

##
# SBayesR
##

if config["testing"] != 'NA':
  n_cores_sbayesr = config.get("ncores", 1)
else:
  n_cores_sbayesr = config.get("ncores", 10)

mem_sbayesr=4000*n_cores_sbayesr

rule prep_pgs_sbayesr_i:
  resources:
    mem_mb=mem_sbayesr
  threads: n_cores_sbayesr
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.download_gctb_ref.output,
    rules.download_gctb_software.output
  output:
    touch(f"{outdir}/reference/target_checks/prep_pgs_sbayesr_i-{{gwas}}.done")
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_sbayesr_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_sbayesr_i-{{gwas}}.log"
  params:
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/sbayesr.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --gctb resources/software/gctb/gctb_2.03beta_Linux/gctb \
      --ld_matrix_chr resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr \
      --robust T \
      --n_cores {threads} \
      --output {outdir}/reference/pgs_score_files/sbayesr/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data resources/data/ref/ref.pop.txt \
      --test {params.testing} > {log} 2>&1"

rule prep_pgs_sbayesr:
  input: expand(f"{outdir}/reference/target_checks/prep_pgs_sbayesr_i-{{gwas}}.done", gwas=gwas_list_df_eur['name'])

##
# lassosum
##

if config["testing"] != 'NA':
  n_cores_lassosum = config.get("ncores", 1)
else:
  n_cores_lassosum = config.get("ncores", 10)

mem_lassosum=10000

rule prep_pgs_lassosum_i:
  resources:
    mem_mb=mem_lassosum
  threads: n_cores_lassosum
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.install_lassosum.output
  output:
    touch(f"{outdir}/reference/target_checks/prep_pgs_lassosum_i-{{gwas}}.done")
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_lassosum_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_lassosum_i-{{gwas}}.log"
  conda:
    "../envs/analysis.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/lassosum.R \
     --ref_plink_chr resources/data/ref/ref.chr \
     --ref_keep resources/data/ref/keep_files/{params.population}.keep \
     --gwas_pop {params.population} \
     --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
     --output {outdir}/reference/pgs_score_files/lassosum/{wildcards.gwas}/ref-{wildcards.gwas} \
      --n_cores {threads} \
     --pop_data resources/data/ref/ref.pop.txt \
     --test {params.testing} > {log} 2>&1"

rule prep_pgs_lassosum:
  input: expand(f"{outdir}/reference/target_checks/prep_pgs_lassosum_i-{{gwas}}.done", gwas=gwas_list_df['name'])

##
# LDpred2
##

if config["testing"] != 'NA':
  n_cores_ldpred2 = config.get("ncores", 5)
else:
  n_cores_ldpred2 = config.get("ncores", 10)

mem_ldpred2=30000

rule prep_pgs_ldpred2_i:
  resources:
    mem_mb=mem_ldpred2,
    time_min=800
  threads: n_cores_ldpred2
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.download_ldpred2_ref.output
  output:
    touch(f"{outdir}/reference/target_checks/prep_pgs_ldpred2_i-{{gwas}}.done")
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_ldpred2_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_ldpred2_i-{{gwas}}.log"
  conda:
    "../envs/analysis.yaml"
  params:
    testing=config["testing"]
  shell:
    "export OPENBLAS_NUM_THREADS=1; \
    Rscript ../Scripts/pgs_methods/ldpred2.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/EUR.keep \
      --ldpred2_ref_dir resources/data/ldpred2_ref \
      --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --n_cores {threads} \
      --output {outdir}/reference/pgs_score_files/ldpred2/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data resources/data/ref/ref.pop.txt \
      --test {params.testing} > {log} 2>&1"

rule prep_pgs_ldpred2:
  input: expand(f"{outdir}/reference/target_checks/prep_pgs_ldpred2_i-{{gwas}}.done", gwas=gwas_list_df_eur['name'])

##
# LDAK MegaPRS
##

if config["testing"] != 'NA':
  n_cores_megaprs = config.get("ncores", 5)
else:
  n_cores_megaprs = config.get("ncores", 10)

mem_megaprs=20000

rule prep_pgs_megaprs_i:
  resources:
    mem_mb=mem_megaprs,
    time_min=800
  threads: n_cores_megaprs
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.download_ldak_highld.output,
    rules.download_ldak.output,
    rules.download_ldak_map.output,
    rules.download_ldak_bld.output
  output:
    touch(f"{outdir}/reference/target_checks/prep_pgs_megaprs_i-{{gwas}}.done")
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_megaprs_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_megaprs_i-{{gwas}}.log"
  conda:
    "../envs/analysis.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/megaprs.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{params.population}.keep \
      --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --ldak resources/software/ldak/ldak5.1.linux \
      --ldak_map resources/data/ldak_map/genetic_map_b37 \
      --ldak_tag resources/data/ldak_bld \
      --ldak_highld resources/data/ldak_highld/highld.txt \
      --n_cores {threads} \
      --output {outdir}/reference/pgs_score_files/megaprs/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data resources/data/ref/ref.pop.txt \
      --test {params.testing} > {log} 2>&1"

rule prep_pgs_megaprs:
  input: expand(f"{outdir}/reference/target_checks/prep_pgs_megaprs_i-{{gwas}}.done", gwas=gwas_list_df['name'])

##
# Process externally created score files
##

# Read in score_list or create empty score_list
if 'score_list' in config and config["score_list"] != 'NA':
  score_list_df = pd.read_table(config["score_list"], sep=r'\s+')
  pgs_methods = config['pgs_methods']
  pgs_methods_all = list(config['pgs_methods'])
  pgs_methods_all.append('external')

  # Check whether score_list paths exist
  check_list_paths(score_list_df)
else:
  score_list_df = pd.DataFrame(columns = ["name", "path", "label"])
  pgs_methods = config['pgs_methods']
  pgs_methods_all = config['pgs_methods']

# Check whether gwas_list paths exist
check_list_paths(score_list_df)

# Download PGS score files for PGSC if path is NA
rule download_pgs_external:
  input:
    rules.download_pgscatalog_utils.output
  output:
    touch(f"{outdir}/reference/pgs_score_files/raw_external/{{score}}/{{score}}_hmPOS_GRCh37.txt.gz")
  params:
    config_file = config["config_file"],
    outdir=config["outdir"],
    path= lambda w: score_list_df.loc[score_list_df['name'] == "{}".format(w.score), 'path'].iloc[0],
    testing = config['testing']
  benchmark:
    f"{outdir}/reference/benchmarks/download_pgs_external-{{score}}.txt"
  log:
    f"{outdir}/reference/logs/download_pgs_external-{{score}}.log"
  conda:
    "../envs/pgscatalog_utils.yaml"
  shell:
    "mkdir -p {outdir}/reference/pgs_score_files/raw_external/{wildcards.score}; \
    download_scorefiles -w -i {wildcards.score} -o {outdir}/reference/pgs_score_files/raw_external/{wildcards.score} -b GRCh37 > {log} 2>&1"

# Create function to return path for score file, depeding on whether the score file was downloaded
def score_path(w):
  if not pd.isna(score_list_df.loc[score_list_df['name'] == w.score, 'path'].iloc[0]):
      return [score_list_df.loc[score_list_df['name'] == w.score, 'path'].iloc[0]]
  else:
      return [outdir + "/reference/pgs_score_files/raw_external/" + w.score + "/" + w.score + "_hmPOS_GRCh37.txt.gz"]

# Harmonise external score files to reference
rule prep_pgs_external_i:
  input:
    lambda w: score_path(w),
    rules.download_default_ref.output,
    rules.install_genoutils.output
  output:
    touch(f"{outdir}/reference/target_checks/prep_pgs_external_i-{{score}}.done")
  params:
    config_file = config["config_file"],
    outdir=config["outdir"],
    score= lambda w: score_path(w),
    testing=config["testing"]
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_external_i-{{score}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_external_i-{{score}}.log"
  conda:
    "../envs/analysis.yaml"
  shell:
    "Rscript ../Scripts/external_score_processor/external_score_processor.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --score {params.score} \
      --output {outdir}/reference/pgs_score_files/external/{wildcards.score}/ref-{wildcards.score} \
      --pop_data resources/data/ref/ref.pop.txt \
      --test {params.testing} > {log} 2>&1"

rule prep_pgs_external:
  input: expand(f"{outdir}/reference/target_checks/prep_pgs_external_i-{{score}}.done", score=score_list_df['name'])

# Create a file listing score files and whether they had sufficient overlap with reference
checkpoint score_reporter:
  input:
    expand(f"{outdir}/reference/target_checks/prep_pgs_external_i-{{score}}.done", score=score_list_df['name'])
  output:
    touch(f"{outdir}/reference/target_checks/score_reporter.done")
  benchmark:
    f"{outdir}/reference/benchmarks/score_reporter.txt"
  log:
    f"{outdir}/reference/logs/score_reporter.log"
  conda:
    "../envs/analysis.yaml"
  params:
    config_file = config["config_file"]
  shell:
    "Rscript ../Scripts/pipeline_misc/score_reporter.R {params.config_file} > {log} 2>&1"

##
# Use a rule to check requested PGS methods have been run for all GWAS
##

pgs_methods_input = list()

if 'ptclump' in pgs_methods_all:
  pgs_methods_input.append(rules.prep_pgs_ptclump.input)
if 'dbslmm' in pgs_methods_all:
  pgs_methods_input.append(rules.prep_pgs_dbslmm.input)
if 'prscs' in pgs_methods_all:
  pgs_methods_input.append(rules.prep_pgs_prscs.input)
if 'sbayesr' in pgs_methods_all:
  pgs_methods_input.append(rules.prep_pgs_sbayesr.input)
if 'lassosum' in pgs_methods_all:
  pgs_methods_input.append(rules.prep_pgs_lassosum.input)
if 'ldpred2' in pgs_methods_all:
  pgs_methods_input.append(rules.prep_pgs_ldpred2.input)
if 'megaprs' in pgs_methods_all:
  pgs_methods_input.append(rules.prep_pgs_megaprs.input)
if 'external' in pgs_methods_all:
  pgs_methods_input.append(rules.prep_pgs_external.input)

rule prep_pgs:
  input:
    pgs_methods_input
