# Create PC score files specific to each population
rule ref_pca_i:
  input:
    ref_input,
    rules.install_genoutils.output,
    f"{resdir}/last_version.txt"
  output:
    f"{resdir}/data/ref/pc_score_files/{{population}}/ref-{{population}}-pcs.EUR.scale"
  conda:
    "../envs/analysis.yaml",
  params:
    testing=config["testing"]
  benchmark:
    f"{resdir}/data/benchmarks/ref_pca_i-{{population}}.txt"
  log:
    f"{resdir}/data/logs/ref_pca_i-{{population}}.log"
  shell:
    "Rscript ../Scripts/ref_pca/ref_pca.R \
      --ref_plink_chr {refdir}/ref.chr \
      --ref_keep {refdir}/keep_files/{wildcards.population}.keep \
      --pop_data {refdir}/ref.pop.txt \
      --output {resdir}/data/ref/pc_score_files/{wildcards.population}/ref-{wildcards.population}-pcs \
      --test {params.testing} > {log} 2>&1"

populations=["AFR","AMR","EAS","EUR","SAS"]

rule ref_pca:
  input: expand(f"{resdir}/data/ref/pc_score_files/{{population}}/ref-{{population}}-pcs.EUR.scale", population=populations)

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
      ref_input,
      rules.install_genoutils.output,
      lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0],
      f"{resdir}/last_version.txt"
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
      refdir=config["refdir"],
      config_file = config["config_file"],
      testing = config['testing'],
      population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
      n= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'n'].iloc[0],
      path= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0]
    shell:
      """
      sumstat_cleaner_script=$(Rscript -e 'cat(system.file("scripts", "sumstat_cleaner.R", package = "GenoUtils"))' 2>&1 | grep -Eo "/.*sumstat_cleaner.R")
      Rscript $sumstat_cleaner_script \
        --sumstats {params.path} \
        --n {params.n} \
        --ref_chr {refdir}/ref.chr \
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
    f"{outdir}/reference/pgs_score_files/ptclump/{{gwas}}/ref-{{gwas}}.score.gz"
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_ptclump_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_ptclump_i-{{gwas}}.log"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    pts= ",".join(map(str, config["ptclump_pts"])),
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/ptclump.R \
      --ref_plink_chr {refdir}/ref.chr \
      --ref_keep {refdir}/keep_files/{params.population}.keep \
      --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --output {outdir}/reference/pgs_score_files/ptclump/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data {refdir}/ref.pop.txt \
      --pTs {params.pts} \
      --test {params.testing} > {log} 2>&1"

rule prep_pgs_ptclump:
  input: expand(f"{outdir}/reference/pgs_score_files/ptclump/{{gwas}}/ref-{{gwas}}.score.gz", gwas=gwas_list_df['name'])

##
# DBSLMM
##

# Create function to specify LD block data from EUR if the GWAS population is CSA, AMR or MID
def set_ld_blocks_pop(population):
  if population in ['CSA', 'AMR', 'MID']:
    return 'EUR'
  else:
    return population

rule prep_pgs_dbslmm_i:
  threads: config['cores_prep_pgs']
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.download_plink.output,
    rules.download_ldscores_panukb.output,
    rules.download_ldsc.output,
    rules.download_hm3_snplist.output,
    rules.download_dbslmm.output,
    rules.download_ld_blocks.output
  output:
    f"{outdir}/reference/pgs_score_files/dbslmm/{{gwas}}/ref-{{gwas}}.score.gz"
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_dbslmm_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_dbslmm_i-{{gwas}}.log"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    ld_block_pop= lambda w: set_ld_blocks_pop(gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0]),
    sampling= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'sampling'].iloc[0],
    prevalence= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'prevalence'].iloc[0],
    h2f= ",".join(map(str, config["dbslmm_h2f"])),
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/dbslmm.R \
      --ref_plink_chr {refdir}/ref.chr \
      --ref_keep {refdir}/keep_files/{params.population}.keep \
      --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --ld_blocks {resdir}/data/ld_blocks/{params.ld_block_pop} \
      --plink {resdir}/software/plink/plink \
      --dbslmm {resdir}/software/dbslmm/software \
      --munge_sumstats {resdir}/software/ldsc/munge_sumstats.py \
      --ldsc {resdir}/software/ldsc/ldsc.py \
      --ld_scores {resdir}/data/ld_scores/UKBB.{params.population}.rsid \
      --hm3_snplist {resdir}/data/hm3_snplist/w_hm3.snplist \
      --hm3_no_mhc T \
      --sample_prev {params.sampling} \
      --pop_prev {params.prevalence} \
      --output {outdir}/reference/pgs_score_files/dbslmm/{wildcards.gwas}/ref-{wildcards.gwas} \
      --n_cores {threads} \
      --pop_data {refdir}/ref.pop.txt \
      --h2f {params.h2f} \
      --test {params.testing} > {log} 2>&1"

rule prep_pgs_dbslmm:
  input: expand(f"{outdir}/reference/pgs_score_files/dbslmm/{{gwas}}/ref-{{gwas}}.score.gz", gwas=gwas_list_df['name'])

##
# PRScs
##
# Note. Threads are set to 1, and phi and chr are run in parallel. Increasing number of threads shows no improvement in speed.

rule prep_pgs_prscs_i:
  resources:
    mem_mb=2000*config['cores_prep_pgs'],
    time_min=800
  threads: config['cores_prep_pgs']
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.download_prscs_software.output,
    lambda w: f"{resdir}/data/prscs_ref/" + prscs_ldref + "/ldblk_" + prscs_ldref + "_" + gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0].lower() + "/ldblk_1kg_chr1.hdf5"
  output:
    f"{outdir}/reference/pgs_score_files/prscs/{{gwas}}/ref-{{gwas}}.score.gz"
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_prscs_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_prscs_i-{{gwas}}.log"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0].lower(),
    phi= ",".join(map(str, config["prscs_phi"])),
    testing=config["testing"]
  shell:
    """
    export MKL_NUM_THREADS=1; \
    export NUMEXPR_NUM_THREADS=1; \
    export OMP_NUM_THREADS=1; \
    export OPENBLAS_NUM_THREADS=1; \
    Rscript ../Scripts/pgs_methods/prscs.R \
    --ref_plink_chr {refdir}/ref.chr \
    --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
    --output {outdir}/reference/pgs_score_files/prscs/{wildcards.gwas}/ref-{wildcards.gwas} \
    --pop_data {refdir}/ref.pop.txt \
    --PRScs_path {resdir}/software/prscs/PRScs.py \
    --PRScs_ref_path {resdir}/data/prscs_ref/{prscs_ldref}/ldblk_{prscs_ldref}_{params.population} \
    --n_cores {threads} \
    --phi_param {params.phi} \
    --test {params.testing} > {log} 2>&1
    """

rule prep_pgs_prscs:
  input: expand(f"{outdir}/reference/pgs_score_files/prscs/{{gwas}}/ref-{{gwas}}.score.gz", gwas=gwas_list_df['name'])

##
# SBayesR
##

rule prep_pgs_sbayesr_i:
  resources:
    mem_mb=4000*config['cores_prep_pgs']
  threads: config['cores_prep_pgs']
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.download_gctb_ref.output,
    rules.download_gctb_software.output
  output:
    f"{outdir}/reference/pgs_score_files/sbayesr/{{gwas}}/ref-{{gwas}}.score.gz"
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
      --ref_plink_chr {refdir}/ref.chr \
      --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --gctb {resdir}/software/gctb/gctb_2.03beta_Linux/gctb \
      --ld_matrix_chr {resdir}/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr \
      --robust T \
      --n_cores {threads} \
      --output {outdir}/reference/pgs_score_files/sbayesr/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data {refdir}/ref.pop.txt \
      --test {params.testing} > {log} 2>&1"

rule prep_pgs_sbayesr:
  input: expand(f"{outdir}/reference/pgs_score_files/sbayesr/{{gwas}}/ref-{{gwas}}.score.gz", gwas=gwas_list_df_eur['name'])

##
# lassosum
##

rule prep_pgs_lassosum_i:
  resources:
    mem_mb=10000
  threads: config['cores_prep_pgs']
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.install_lassosum.output
  output:
    f"{outdir}/reference/pgs_score_files/lassosum/{{gwas}}/ref-{{gwas}}.score.gz"
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
     --ref_plink_chr {refdir}/ref.chr \
     --ref_keep {refdir}/keep_files/{params.population}.keep \
     --gwas_pop {params.population} \
     --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
     --output {outdir}/reference/pgs_score_files/lassosum/{wildcards.gwas}/ref-{wildcards.gwas} \
     --n_cores {threads} \
     --pop_data {refdir}/ref.pop.txt \
     --test {params.testing} > {log} 2>&1"

rule prep_pgs_lassosum:
  input: expand(f"{outdir}/reference/pgs_score_files/lassosum/{{gwas}}/ref-{{gwas}}.score.gz", gwas=gwas_list_df['name'])

##
# LDpred2
##

rule prep_pgs_ldpred2_i:
  resources:
    mem_mb=30000,
    time_min=800
  threads: config['cores_prep_pgs']
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.download_ldpred2_ref.output
  output:
    f"{outdir}/reference/pgs_score_files/ldpred2/{{gwas}}/ref-{{gwas}}.score.gz"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_ldpred2_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_ldpred2_i-{{gwas}}.log"
  conda:
    "../envs/analysis.yaml"
  params:
    model=",".join(map(str, config["ldpred2_model"])),
    inference=",".join(map(str, config["ldpred2_inference"])),
    sampling= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'sampling'].iloc[0],
    binary=lambda w: 'T' if not pd.isna(gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'sampling'].iloc[0]) else 'F',
    testing=config["testing"]
  shell:
    "export OPENBLAS_NUM_THREADS=1; \
    Rscript ../Scripts/pgs_methods/ldpred2.R \
      --ref_plink_chr {refdir}/ref.chr \
      --ref_keep {refdir}/keep_files/EUR.keep \
      --ldpred2_ref_dir {resdir}/data/ldpred2_ref \
      --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --n_cores {threads} \
      --output {outdir}/reference/pgs_score_files/ldpred2/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data {refdir}/ref.pop.txt \
      --model {params.model} \
      --binary {params.binary} \
      --inference {params.inference} \
      --sample_prev {params.sampling} \
      --test {params.testing} > {log} 2>&1"

rule prep_pgs_ldpred2:
  input: expand(f"{outdir}/reference/pgs_score_files/ldpred2/{{gwas}}/ref-{{gwas}}.score.gz", gwas=gwas_list_df_eur['name'])

##
# LDAK MegaPRS
##

rule prep_pgs_megaprs_i:
  resources:
    mem_mb=20000,
    time_min=800
  threads: config['cores_prep_pgs']
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.download_ldak_highld.output,
    rules.download_ldak.output,
    rules.download_ldak_map.output,
    rules.download_ldak_bld.output
  output:
    f"{outdir}/reference/pgs_score_files/megaprs/{{gwas}}/ref-{{gwas}}.score.gz"
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
      --ref_plink_chr {refdir}/ref.chr \
      --ref_keep {refdir}/keep_files/{params.population}.keep \
      --sumstats {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --ldak {resdir}/software/ldak/ldak5.1.linux \
      --ldak_map {resdir}/data/ldak_map/genetic_map_b37 \
      --ldak_tag {resdir}/data/ldak_bld \
      --ldak_highld {resdir}/data/ldak_highld/highld.txt \
      --n_cores {threads} \
      --output {outdir}/reference/pgs_score_files/megaprs/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data {refdir}/ref.pop.txt \
      --test {params.testing} > {log} 2>&1"

rule prep_pgs_megaprs:
  input: expand(f"{outdir}/reference/pgs_score_files/megaprs/{{gwas}}/ref-{{gwas}}.score.gz", gwas=gwas_list_df['name'])

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
    rules.download_pgscatalog_utils.output,
    f"{resdir}/last_version.txt"
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
    pgscatalog-download --pgs {wildcards.score} -o {outdir}/reference/pgs_score_files/raw_external/{wildcards.score} --build GRCh37 > {log} 2>&1"

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
    ref_input,
    rules.install_genoutils.output
  output:
    touch(f"{outdir}/reference/target_checks/prep_pgs_external_i-{{score}}.done")
  params:
    config_file = config["config_file"],
    outdir=config["outdir"],
    refdir=config["refdir"],
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
      --ref_plink_chr {refdir}/ref.chr \
      --score {params.score} \
      --output {outdir}/reference/pgs_score_files/external/{wildcards.score}/ref-{wildcards.score} \
      --pop_data {refdir}/ref.pop.txt \
      --test {params.testing} > {log} 2>&1"

rule prep_pgs_external:
  input: expand(f"{outdir}/reference/pgs_score_files/external/{{score}}/ref-{{score}}.log", score=score_list_df['name'])

# Create a file listing score files and whether they had sufficient overlap with reference
checkpoint score_reporter:
  input:
    expand(f"{outdir}/reference/target_checks/prep_pgs_external_i-{{score}}.done", score=score_list_df['name'])
  output:
    f"{outdir}/reference/pgs_score_files/external/score_report.txt"
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

###########
# Multi-ancestry methods
###########

# Read in the gwas_groups or make an empty version
if 'gwas_groups' in config and config["gwas_groups"] != 'NA':
  gwas_groups_df = pd.read_table(config["gwas_groups"], sep=r'\s+')
else:
  gwas_groups_df = pd.DataFrame(columns = ["name", "gwas", "label"])

# Function to get the list of GWAS names for a given group
def get_gwas_names(gwas_group):
    gwas_names_str = gwas_groups_df[gwas_groups_df['name'] == gwas_group]['gwas'].iloc[0]
    return gwas_names_str.split(',')

# Function to generate comma-separated list of populations for each name
def get_populations(gwas_group):
    gwas_names = get_gwas_names(gwas_group)
    sumstats_populations = []
    for gwas in gwas_names:
        gwas_info = gwas_list_df[gwas_list_df['name'] == gwas].iloc[0]
        sumstats_populations.append(gwas_info['population'])
    return sumstats_populations

####
# PRS-CSx
####

# Note. Threads are set to 1, and phi and chr are run in parallel. Increasing number of threads shows no improvement in speed.

rule prep_pgs_prscsx_i:
  resources:
    mem_mb=2000*config['cores_prep_pgs'],
    time_min=800
  threads: config['cores_prep_pgs']
  input:
    rules.download_prscsx_software.output,
    lambda w: expand(f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz", gwas=get_gwas_names(w.gwas_group)),
    lambda w: expand(f"{resdir}/data/prscs_ref/{prscs_ldref}/ldblk_{prscs_ldref}_{{population}}/ldblk_{prscs_ldref}_chr1.hdf5", population=[pop.lower() for pop in get_populations(w.gwas_group)]),
    f"{resdir}/data/prscs_ref/{prscs_ldref}/snpinfo_mult_{prscs_ldref}_hm3"
  output:
    f"{outdir}/reference/pgs_score_files/prscsx/{{gwas_group}}/ref-{{gwas_group}}.score.gz"
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_prscsx_i-{{gwas_group}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_prscsx_i-{{gwas_group}}.log"
  params:
    sumstats= lambda w: ",".join(expand(f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz", gwas=get_gwas_names(w.gwas_group))),
    populations= lambda w: ",".join(get_populations(w.gwas_group)),
    phi= ",".join(map(str, config["prscs_phi"])),
    testing=config["testing"]
  shell:
    """
    export MKL_NUM_THREADS=1; \
    export NUMEXPR_NUM_THREADS=1; \
    export OMP_NUM_THREADS=1; \
    export OPENBLAS_NUM_THREADS=1; \
    Rscript ../Scripts/pgs_methods/prscsx.R \
      --ref_plink_chr {refdir}/ref.chr \
      --sumstats {params.sumstats} \
      --populations {params.populations} \
      --prscsx_ref_path {resdir}/data/prscs_ref/{prscs_ldref} \
      --phi_param {params.phi} \
      --pop_data {refdir}/ref.pop.txt \
      --prscsx_path {resdir}/software/prscsx/PRScsx.py \
      --output {outdir}/reference/pgs_score_files/prscsx/{wildcards.gwas_group}/ref-{wildcards.gwas_group} \
      --test {params.testing} \
      --n_cores {threads} > {log} 2>&1
    """

rule prep_pgs_prscsx:
  input: expand(f"{outdir}/reference/pgs_score_files/prscsx/{{gwas_group}}/ref-{{gwas_group}}.score.gz", gwas_group=gwas_groups_df['name'])

####
# X-WING
####

rule prep_pgs_xwing_i:
  resources:
    mem_mb=2000*config['cores_prep_pgs'],
    time_min=800
  threads: config['cores_prep_pgs']
  input:
    rules.download_xwing_software.output,
    rules.download_leopard_panther_snp_data.output,
    lambda w: expand(f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz", gwas=get_gwas_names(w.gwas_group)),
    lambda w: expand(f"{resdir}/data/prscs_ref/{prscs_ldref}/ldblk_{prscs_ldref}_{{population}}/ldblk_{prscs_ldref}_chr1.hdf5", population=[pop.lower() for pop in get_populations(w.gwas_group)]),
    f"{resdir}/data/prscs_ref/{prscs_ldref}/snpinfo_mult_{prscs_ldref}_hm3",
    lambda w: expand(f"{resdir}/data/PANTHER_LEOPARD_1kg_ref/ldblk_1kg_{{population}}/ldblk_1kg_chr13.hdf5", population=[pop.lower() for pop in get_populations(w.gwas_group)]),
    lambda w: expand(f"{resdir}/data/LEOPARD_1kg_ref/{{population}}/{{population}}_part1.bed", population=get_populations(w.gwas_group)),
    lambda w: expand(f"{resdir}/data/LOGODetect_1kg_ref/{{population}}/1000G_{{population}}_QC.bim", population=get_populations(w.gwas_group))
  output:
    f"{outdir}/reference/pgs_score_files/xwing/{{gwas_group}}/ref-{{gwas_group}}.score.gz"
  conda:
    "../envs/xwing.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_xwing_i-{{gwas_group}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_xwing_i-{{gwas_group}}.log"
  params:
    sumstats= lambda w: ",".join(expand(f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz", gwas=get_gwas_names(w.gwas_group))),
    populations= lambda w: ",".join(get_populations(w.gwas_group)),
    testing=config["testing"]
  shell:
    """
    export MKL_NUM_THREADS=1; \
    export NUMEXPR_NUM_THREADS=1; \
    export OMP_NUM_THREADS=1; \
    export OPENBLAS_NUM_THREADS=1; \
    Rscript ../Scripts/pgs_methods/xwing.R \
      --ref_plink_chr {refdir}/ref.chr \
      --sumstats {params.sumstats} \
      --populations {params.populations} \
      --logodetect_ref {resdir}/data/LOGODetect_1kg_ref \
      --panther_ref {resdir}/data/prscs_ref/{prscs_ldref} \
      --leopard_ref {resdir}/data/LEOPARD_1kg_ref \
      --panther_leopard_ref {resdir}/data/PANTHER_LEOPARD_1kg_ref \
      --xwing_repo {resdir}/software/xwing \
      --pop_data {refdir}/ref.pop.txt \
      --output {outdir}/reference/pgs_score_files/xwing/{wildcards.gwas_group}/ref-{wildcards.gwas_group} \
      --test {params.testing} \
      --n_cores {threads} > {log} 2>&1
    """

rule prep_pgs_xwing:
  input: expand(f"{outdir}/reference/pgs_score_files/xwing/{{gwas_group}}/ref-{{gwas_group}}.score.gz", gwas_group=gwas_groups_df['name'])

####
# TL-PRS
####

rule prep_pgs_tlprs_i:
  resources:
    mem_mb=10000,
    time_min=800
  threads: config['cores_prep_pgs']
  input:
    rules.install_tlprs.output,
    lambda w: expand(f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz", gwas=get_gwas_names(w.gwas_group)),
    lambda w: expand(f"{outdir}/reference/pgs_score_files/{{method}}/{{gwas}}/ref-{{gwas}}.score.gz", gwas=get_gwas_names(w.gwas_group), method=w.method)
  output:
    f"{outdir}/reference/pgs_score_files/tlprs_{{method}}/{{gwas_group}}/ref-{{gwas_group}}.score.gz"
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_tlprs_i-{{gwas_group}}-{{method}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_tlprs_i-{{gwas_group}}-{{method}}.log"
  params:
    sumstats= lambda w: ",".join(expand(f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz", gwas=get_gwas_names(w.gwas_group))),
    scores= lambda w: ",".join(expand(f"{outdir}/reference/pgs_score_files/{{method}}/{{gwas}}/ref-{{gwas}}.score.gz", gwas=get_gwas_names(w.gwas_group), method=w.method)),
    populations= lambda w: ",".join(get_populations(w.gwas_group)),
    testing=config["testing"]
  shell:
    """
    Rscript ../Scripts/pgs_methods/tlprs.R \
      --ref_plink_chr {refdir}/ref.chr \
      --sumstats {params.sumstats} \
      --scores {params.scores} \
      --populations {params.populations} \
      --pop_data {refdir}/ref.pop.txt \
      --ref_keep_dir {refdir}/keep_files \
      --output {outdir}/reference/pgs_score_files/tlprs_{wildcards.method}/{wildcards.gwas_group}/ref-{wildcards.gwas_group} \
      --test {params.testing} \
      --n_cores {threads} > {log} 2>&1
    """

rule prep_pgs_tlprs:
  input: expand(f"{outdir}/reference/pgs_score_files/tlprs_{{method}}/{{gwas_group}}/ref-{{gwas_group}}.score.gz", gwas_group=gwas_groups_df['name'], method=config["tlprs_methods"])

####
# BridgePRS
####

rule prep_pgs_bridgeprs_i:
  resources:
    mem_mb=2000*config['cores_prep_pgs'],
    time_min=800
  threads: config['cores_prep_pgs']
  input:
    rules.download_bridgeprs_software.output,
    lambda w: expand(f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz", gwas=get_gwas_names(w.gwas_group))
  output:
    f"{outdir}/reference/pgs_score_files/bridgeprs/{{gwas_group}}/ref-{{gwas_group}}.score.gz"
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_bridgeprs_i-{{gwas_group}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_bridgeprs_i-{{gwas_group}}.log"
  params:
    sumstats= lambda w: ",".join(expand(f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz", gwas=get_gwas_names(w.gwas_group))),
    populations= lambda w: ",".join(get_populations(w.gwas_group)),
    testing=config["testing"]
  shell:
    """
    Rscript ../Scripts/pgs_methods/bridgeprs.R \
      --ref_plink_chr {refdir}/ref.chr \
      --sumstats {params.sumstats} \
      --populations {params.populations} \
      --pop_data {refdir}/ref.pop.txt \
      --output {outdir}/reference/pgs_score_files/bridgeprs/{wildcards.gwas_group}/ref-{wildcards.gwas_group} \
      --test {params.testing} \
      --bridgeprs_repo {resdir}/software/bridgeprs \
      --n_cores {threads} > {log} 2>&1
    """

rule prep_pgs_bridgeprs:
  input: expand(f"{outdir}/reference/pgs_score_files/bridgeprs/{{gwas_group}}/ref-{{gwas_group}}.score.gz", gwas_group=gwas_groups_df['name'])

###############################################

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
  pgs_methods_input.append(rules.score_reporter.output)
if 'prscsx' in pgs_methods_all:
  pgs_methods_input.append(rules.prep_pgs_prscsx.input)
if 'xwing' in pgs_methods_all:
  pgs_methods_input.append(rules.prep_pgs_xwing.input)
if 'tlprs' in pgs_methods_all:
  pgs_methods_input.append(rules.prep_pgs_tlprs.input)
if 'bridgeprs' in pgs_methods_all:
  pgs_methods_input.append(rules.prep_pgs_bridgeprs.input)

rule prep_pgs:
  input:
    pgs_methods_input
