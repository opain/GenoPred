# Create PC score files specific to each population
rule ref_pca:
  input:
    rules.get_dependencies.output,
    "../Scripts/ref_pca/ref_pca.R",
    "../Scripts/functions/misc.R"
  output:
    "resources/data/ref/pc_score_files/{population}/ref.{population}.pcs.EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  shell:
    "Rscript ../Scripts/ref_pca/ref_pca.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{wildcards.population}.keep \
      --pop_data resources/data/ref/ref.pop.txt \
      --output resources/data/ref/pc_score_files/{wildcards.population}/ref-{wildcards.population}-pcs"

populations=["AFR","AMR","EAS","EUR","SAS"]

rule run_ref_pca:
  input: expand("resources/data/ref/pc_score_files/{population}/ref-{population}-pcs.EUR.scale", population=populations)
  
##
# QC and format GWAS summary statistics
##

if 'gwas_list' in config:
  gwas_list_df = pd.read_table(config["gwas_list"], sep=r'\s+')
else:
  gwas_list_df = pd.DataFrame(columns = ["name", "path", "population", "n", "sampling", "prevalence", "mean", "sd", "label"])

gwas_list_df_eur = gwas_list_df.loc[gwas_list_df['population'] == 'EUR']

if 'gwas_list' in config:
  rule sumstat_prep:
    input:
      config['gwas_list'],
      config['config_file'],
      rules.get_dependencies.output
    output:
      "{outdir}/resources/data/gwas_sumstat/{gwas}/{gwas}-cleaned.gz"
    conda:
      "../envs/GenoPredPipe.yaml"
    params:
      population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
      path= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0]
    shell:
      """
      sumstat_cleaner_script=$(Rscript -e 'cat(system.file("scripts", "sumstat_cleaner.R", package = "GenoUtils"))')
      Rscript $sumstat_cleaner_script \
        --sumstats {params.path} \
        --ref_chr resources/data/ref/ref.chr \
        --population {params.population} \
        --output {outdir}/resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned
      """
          
  rule run_sumstat_prep:
    input: expand("{outdir}/resources/data/gwas_sumstat/{gwas}/{gwas}-cleaned.gz", gwas=gwas_list_df['name'], outdir=outdir)

##
# pT+clump (sparse, nested)
##

rule prep_pgs_ptclump:
  input:
    "{outdir}/resources/data/gwas_sumstat/{gwas}/{gwas}-cleaned.gz",
    "../Scripts/pgs_methods/ptclump.R",
    "../Scripts/functions/misc.R"
  output:
    "{outdir}/resources/data/ref/pgs_score_files/ptclump/{gwas}/ref-{gwas}-EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/ptclump.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{params.population}.keep \
      --sumstats {outdir}/resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --output {outdir}/resources/data/ref/pgs_score_files/ptclump/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data resources/data/ref/ref.pop.txt \
      --test {params.testing}"
  
rule run_prep_pgs_ptclump:
  input: expand("{outdir}/resources/data/ref/pgs_score_files/ptclump/{gwas}/ref-{gwas}-EUR.scale", gwas=gwas_list_df['name'], outdir=outdir)

##
# DBSLMM
##

rule prep_pgs_dbslmm:
  input:
    "{outdir}/resources/data/gwas_sumstat/{gwas}/{gwas}-cleaned.gz",
    rules.get_dependencies.output,
    "../Scripts/pgs_methods/dbslmm.R",
    "../Scripts/functions/misc.R"
  output:
    "{outdir}/resources/data/ref/pgs_score_files/dbslmm/{gwas}/ref-{gwas}-EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    sampling= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'sampling'].iloc[0],
    prevalence= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'prevalence'].iloc[0],
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/dbslmm.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{params.population}.keep \
      --sumstats {outdir}/resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --ld_blocks resources/data/ld_blocks/{params.population} \
      --plink resources/software/plink/plink \
      --dbslmm resources/software/dbslmm/software \
      --munge_sumstats resources/software/ldsc/munge_sumstats.py \
      --ldsc resources/software/ldsc/ldsc.py \
      --ldsc_ref resources/data/ldsc_ref/eur_w_ld_chr \
      --hm3_snplist resources/data/hm3_snplist/w_hm3.snplist \
      --sample_prev {params.sampling} \
      --pop_prev {params.prevalence} \
      --output {outdir}/resources/data/ref/pgs_score_files/dbslmm/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data resources/data/ref/ref.pop.txt \
      --test {params.testing}"

rule run_prep_pgs_dbslmm:
  input: expand("{outdir}/resources/data/ref/pgs_score_files/dbslmm/{gwas}/ref-{gwas}-EUR.scale", gwas=gwas_list_df_eur['name'], outdir=outdir)

##
# PRScs
##

if config["testing"] != 'NA':
  n_cores_prscs=min(5, multiprocessing.cpu_count())
  mem_prscs=40000
else:
  n_cores_prscs=min(10, multiprocessing.cpu_count())
  mem_prscs=80000

rule prep_pgs_prscs:
  threads:n_cores_prscs
  resources: 
    mem_mb=mem_prscs,
    cpus=n_cores_prscs,
    time_min=800
  input:
    "{outdir}/resources/data/gwas_sumstat/{gwas}/{gwas}-cleaned.gz",
    rules.get_dependencies.output,
    "../Scripts/pgs_methods/prscs.R",
    "../Scripts/functions/misc.R",
    rules.download_prscs_software.output,
    rules.download_prscs_ref_1kg_eur.output
  output:
    "{outdir}/resources/data/ref/pgs_score_files/prscs/{gwas}/ref-{gwas}-EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    testing=config["testing"]
  shell:
    "export MKL_NUM_THREADS=1; \
     export NUMEXPR_NUM_THREADS=1; \
     export OMP_NUM_THREADS=1; \
     Rscript ../Scripts/pgs_methods/prscs.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --sumstats {outdir}/resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --output {outdir}/resources/data/ref/pgs_score_files/prscs/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data resources/data/ref/ref.pop.txt \
      --PRScs_path resources/software/prscs/PRScs.py \
      --PRScs_ref_path resources/data/prscs_ref/ldblk_1kg_eur \
      --n_cores {n_cores_prscs} \
      --phi_param 1e-6,1e-4,1e-2,1,auto \
      --test {params.testing}"

rule run_prep_pgs_prscs:
  input: expand("{outdir}/resources/data/ref/pgs_score_files/prscs/{gwas}/ref-{gwas}-EUR.scale", gwas=gwas_list_df_eur['name'], outdir=outdir)

##
# SBayesR
##

if config["testing"] != 'NA':
  n_cores_sbayesr=min(1, multiprocessing.cpu_count())
  mem_sbayesr=10000
else:
  n_cores_sbayesr=min(10, multiprocessing.cpu_count())
  mem_sbayesr=80000

rule prep_pgs_sbayesr:
  threads:n_cores_sbayesr
  resources: 
    mem_mb=mem_sbayesr,
    cpus=n_cores_sbayesr
  input:
    "{outdir}/resources/data/gwas_sumstat/{gwas}/{gwas}-cleaned.gz",
    rules.get_dependencies.output,
    "../Scripts/pgs_methods/sbayesr.R",
    "../Scripts/functions/misc.R",
    rules.download_gctb_ref.output,
    rules.download_gctb_software.output
  output:
    "{outdir}/resources/data/ref/pgs_score_files/sbayesr/{gwas}/ref-{gwas}-EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/sbayesr.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --sumstats {outdir}/resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --gctb resources/software/gctb/gctb_2.03beta_Linux/gctb \
      --ld_matrix_chr resources/data/gctb_ref/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_v3_50k_chr \
      --robust T \
      --n_cores {n_cores_sbayesr} \
      --output {outdir}/resources/data/ref/pgs_score_files/sbayesr/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data resources/data/ref/ref.pop.txt \
      --test {params.testing}"

rule run_prep_pgs_sbayesr:
  input: expand("{outdir}/resources/data/ref/pgs_score_files/sbayesr/{gwas}/ref-{gwas}-EUR.scale", gwas=gwas_list_df_eur['name'], outdir=outdir)

##
# lassosum
##

rule prep_pgs_lassosum:
  input:
    "{outdir}/resources/data/gwas_sumstat/{gwas}/{gwas}-cleaned.gz",
    rules.get_dependencies.output,
    "../Scripts/pgs_methods/lassosum.R",
    "../Scripts/functions/misc.R"
  output:
    "{outdir}/resources/data/ref/pgs_score_files/lassosum/{gwas}/ref-{gwas}-EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/lassosum.R \
     --ref_plink_chr resources/data/ref/ref.chr \
     --ref_keep resources/data/ref/keep_files/{params.population}.keep \
     --gwas_pop {params.population} \
     --sumstats {outdir}/resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
     --output {outdir}/resources/data/ref/pgs_score_files/lassosum/{wildcards.gwas}/ref-{wildcards.gwas} \
     --pop_data resources/data/ref/ref.pop.txt \
     --test {params.testing}"
    
rule run_prep_pgs_lassosum:
  input: expand("{outdir}/resources/data/ref/pgs_score_files/lassosum/{gwas}/ref-{gwas}-EUR.scale", gwas=gwas_list_df['name'], outdir=outdir)

##
# LDpred2
##

if config["testing"] != 'NA':
  n_cores_ldpred2=min(5, multiprocessing.cpu_count())
  mem_ldpred2=40000
else:
  n_cores_ldpred2=min(10, multiprocessing.cpu_count())
  mem_ldpred2=80000

rule prep_pgs_ldpred2:
  threads:n_cores_ldpred2
  resources: 
    mem_mb=mem_ldpred2,
    cpus=n_cores_ldpred2,
    time_min=800
  input:
    rules.get_dependencies.output,
    "{outdir}/resources/data/gwas_sumstat/{gwas}/{gwas}-cleaned.gz",
    "../Scripts/pgs_methods/ldpred2.R",
    "../Scripts/functions/misc.R",
    rules.download_ldpred2_ref.output
  output:
    "{outdir}/resources/data/ref/pgs_score_files/ldpred2/{gwas}/ref-{gwas}-EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    testing=config["testing"]
  shell:
    "export OPENBLAS_NUM_THREADS=1; \
    Rscript ../Scripts/pgs_methods/ldpred2.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/EUR.keep \
      --ldpred2_ref_dir resources/data/ldpred2_ref \
      --sumstats {outdir}/resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --n_cores {n_cores_ldpred2} \
      --output {outdir}/resources/data/ref/pgs_score_files/ldpred2/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data resources/data/ref/ref.pop.txt \
      --test {params.testing}"
    
rule run_prep_pgs_ldpred2:
  input: expand("{outdir}/resources/data/ref/pgs_score_files/ldpred2/{gwas}/ref-{gwas}-EUR.scale", gwas=gwas_list_df_eur['name'], outdir=outdir)

##
# LDAK MegaPRS
##

if config["testing"] != 'NA':
  n_cores_megaprs=min(5, multiprocessing.cpu_count())
  mem_megaprs=40000
else:
  n_cores_megaprs=min(10, multiprocessing.cpu_count())
  mem_megaprs=80000

rule prep_pgs_megaprs:
  threads:n_cores_megaprs
  resources: 
    mem_mb=mem_megaprs,
    cpus=n_cores_megaprs,
    time_min=800
  input:
    rules.get_dependencies.output,
    "{outdir}/resources/data/gwas_sumstat/{gwas}/{gwas}-cleaned.gz",
    "../Scripts/pgs_methods/megaprs.R",
    "../Scripts/functions/misc.R",
    rules.download_ldak_highld.output,
    rules.download_ldak.output,
    rules.download_ldak_bld.output
  output:
    "{outdir}/resources/data/ref/pgs_score_files/megaprs/{gwas}/ref-{gwas}-EUR.scale"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    testing=config["testing"]
  shell:
    "Rscript ../Scripts/pgs_methods/megaprs.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{params.population}.keep \
      --sumstats {outdir}/resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --ldak resources/software/ldak/ldak5.1.linux \
      --ldak_map resources/data/ldak_map/genetic_map_b37 \
      --ldak_tag resources/data/ldak_bld \
      --ldak_highld resources/data/ldak_highld/highld.txt \
      --n_cores {n_cores_megaprs} \
      --output {outdir}/resources/data/ref/pgs_score_files/megaprs/{wildcards.gwas}/ref-{wildcards.gwas} \
      --pop_data resources/data/ref/ref.pop.txt \
      --test {params.testing}"
    
rule run_prep_pgs_megaprs:
  input: expand("{outdir}/resources/data/ref/pgs_score_files/megaprs/{gwas}/ref-{gwas}-EUR.scale", gwas=gwas_list_df['name'], outdir=outdir)

##
# Process externally created score files
##

if 'score_list' in config and config["score_list"] != 'NA':
  score_list_df = pd.read_table(config["score_list"], sep=r'\s+')
  pgs_methods = config['pgs_methods'] 
  pgs_methods.append('external')
else:
  score_list_df = pd.DataFrame(columns = ["name", "path", "label"])
  pgs_methods = config['pgs_methods']

if 'score_list' in config:
  rule prep_pgs_external:
    input:
      config['score_list'],
      rules.get_dependencies.output,
      lambda w: score_list_df.loc[score_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0],
      "../Scripts/external_score_processor/external_score_processor.R",
      "../Scripts/functions/misc.R"
    output:
      "{outdir}/resources/data/ref/pgs_score_files/external/{gwas}/ref-{gwas}-EUR.scale"
    params:
      score= lambda w: score_list_df.loc[score_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0]
    conda:
      "../envs/GenoPredPipe.yaml"
    shell:
      "Rscript ../Scripts/external_score_processor/external_score_processor.R \
        --ref_plink_chr resources/data/ref/ref.chr \
        --score {params.score} \
        --plink2 plink2 \
        --output {outdir}/resources/data/ref/pgs_score_files/external/{wildcards.gwas}/ref-{wildcards.gwas} \
        --ref_pop_scale resources/data/ref/ref.keep.list"
      
  rule run_prep_pgs_external:
    input: expand("{outdir}/resources/data/ref/pgs_score_files/external/{gwas}/ref-{gwas}-EUR.scale", gwas=score_list_df['name'], outdir=outdir)

##
# Use a rule to check requested PGS methods have been run for all GWAS
##

pgs_methods_input = list()

if 'ptclump' in pgs_methods:
  pgs_methods_input.append(rules.run_prep_pgs_ptclump.input)
if 'dbslmm' in pgs_methods:
  pgs_methods_input.append(rules.run_prep_pgs_dbslmm.input)
if 'prscs' in pgs_methods:
  pgs_methods_input.append(rules.run_prep_pgs_prscs.input)
if 'sbayesr' in pgs_methods:
  pgs_methods_input.append(rules.run_prep_pgs_sbayesr.input)
if 'lassosum' in pgs_methods:
  pgs_methods_input.append(rules.run_prep_pgs_lassosum.input)
if 'ldpred2' in pgs_methods:
  pgs_methods_input.append(rules.run_prep_pgs_ldpred2.input)
if 'megaprs' in pgs_methods:
  pgs_methods_input.append(rules.run_prep_pgs_megaprs.input)
if 'external' in pgs_methods:
  pgs_methods_input.append(rules.run_prep_pgs_external.input)

rule pgs_methods_complete:
  input:
    pgs_methods_input
  
##
# Estimate R2/AUC of PRS using lassosum pseudovalidate
##

rule pseudovalidate_prs:
  input:
    "{outdir}/resources/data/gwas_sumstat/{gwas}/{gwas}-cleaned.gz",
    rules.get_dependencies.output,
    "../Scripts/lassosum_pseudovalidate/lassosum_pseudovalidate.R"
  output:
    "{outdir}/resources/data/ref/prs_pseudoval/{gwas}/lassosum_pseudo-{gwas}-pseudovalidate.png"
  conda:
    "../envs/GenoPredPipe.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
  shell:
    "Rscript ../Scripts/lassosum_pseudovalidate/lassosum_pseudovalidate.R \
      --ref_plink_chr resources/data/ref/ref.chr \
      --ref_keep resources/data/ref/keep_files/{params.population}.keep \
      --sumstats {outdir}/resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz \
      --prune_mhc T \
      --output {outdir}/resources/data/ref/prs_pseudoval/{wildcards.gwas}/lassosum_pseudo-{wildcards.gwas} \
      --plink plink \
      --n_cores 1"

rule run_pseudovalidate_prs:
  input: expand("{outdir}/resources/data/ref/prs_pseudoval/{gwas}/lassosum_pseudo-{gwas}-pseudovalidate.png", gwas=gwas_list_df['name'], outdir=outdir)

