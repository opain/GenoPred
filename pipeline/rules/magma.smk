# This code has been adapted from GenoDisc

####
# Download MAGMA
####

rule download_magma:
  output:
    f"{resdir}/software/magma/magma"
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/download_magma.txt"
  log:
    f"{outdir}/reference/logs/download_magma.log"
  shell:
    "wget --no-check-certificate -O {resdir}/software/magma.zip https://vu.data.surfsara.nl/index.php/s/zkKbNeNOZAhFXZB/download; \
     unzip resources/software/magma.zip -d {resdir}/software/magma; \
     rm {resdir}/software/magma.zip > {log} 2>&1"

####
# Download MAGMA gene locations
####

rule download_magma_gene_loc:
  output:
    f"{resdir}/data/magma/NCBI37.3.gene.loc"
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/download_magma_gene_loc.txt"
  log:
    f"{outdir}/reference/logs/download_magma_gene_loc.log"
  shell:
    "wget --no-check-certificate -O {resdir}/data/magma.zip https://vu.data.surfsara.nl/index.php/s/Pj2orwuF2JYyKxq/download; \
     unzip {resdir}/data/magma.zip -d {resdir}/data/magma; \
     rm {resdir}/data/magma.zip > {log} 2>&1"

####
# Download MAGMA reference
####

rule download_magma_ref:
  output:
    f"{resdir}/data/magma_ref/g1000_eur.bed"
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/download_magma_ref.txt"
  log:
    f"{outdir}/reference/logs/download_magma_ref.log"
  shell:
    "wget --no-check-certificate -O {resdir}/data/magma.zip https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download; \
     unzip {resdir}/data/magma.zip -d {resdir}/data/magma_ref; \
     rm {resdir}/data/magma.zip > {log} 2>&1"

####
# Create MAGMA annotation file
####

rule magma_annot:
  input:
    rules.download_magma.output,
    rules.download_magma_gene_loc.output,
    rules.download_magma_ref.output
  output:
    f"{resdir}/data/magma/NCBI37.3.genes.annot"
  conda: 
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/magma_annot.txt"
  log:
    f"{outdir}/reference/logs/magma_annot.log"
  shell:
    "{resdir}/software/magma/magma \
      --annotate window=35,10 \
      --snp-loc {resdir}/data/magma_ref/g1000_eur.bim \
      --gene-loc {resdir}/data/magma/NCBI37.3.gene.loc \
      --out {resdir}/data/magma/NCBI37.3 > {log} 2>&1"

####
# MAGMA
####

# Gene level association analysis
rule magma_gene_level:
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    rules.magma_annot.output
  output:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/magma/magma_gene_level.genes.raw"
  conda: 
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/magma_gene_level-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/magma_gene_level-{{gwas}}.log"
  shell:
    "gzip -f -d -c {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned.gz > {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned; \
     mkdir -p {outdir}/reference/gwas_sumstat/{wildcards.gwas}/magma/; \
     resources/software/magma/magma \
      --bfile {resdir}/data/magma_ref/g1000_eur \
      --pval {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned use=SNP,P ncol=N \
      --gene-annot {resdir}/data/magma/NCBI37.3.genes.annot \
      --out {outdir}/reference/gwas_sumstat/{wildcards.gwas}/magma/magma_gene_level; \
     rm {outdir}/reference/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}-cleaned > {log} 2>&1"

######
# Gene set enrichment analysis
######

rule magma_gene_set_level:
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/magma/magma_gene_level.genes.raw"
  output:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/magma/magma_set_level.gsa.out"
  conda: 
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/magma_gene_set_level-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/magma_gene_set_level-{{gwas}}.log"
  params:
    gene_sets=config['gene_sets']
  shell:
    "resources/software/magma/magma \
      --gene-results {outdir}/reference/gwas_sumstat/{wildcards.gwas}/magma/magma_gene_level.genes.raw \
      --set-annot {params.gene_sets} \
      --model direction-sets=greater \
      --out {outdir}/reference/gwas_sumstat/{wildcards.gwas}/magma/magma_set_level > {log} 2>&1"

# Run conditional analysis
rule magma_set_conditional:
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/magma/magma_set_level.gsa.out"
  output:
    touch(f"{outdir}/reference/gwas_sumstat/{{gwas}}/magma/magma_set_conditional.done")
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/magma_set_conditional-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/magma_set_conditional-{{gwas}}.log"
  params:
    config_file= config['config_file']
  shell:
    "Rscript ../Scripts/magma/magma_set_conditional.R \
      --config {params.config_file} \
      --gwas {wildcards.gwas} > {log} 2>&1"

# Create SNP-lists for enriched gene sets
rule create_set_snplists:
  input:
    f"{outdir}/reference/gwas_sumstat/{{gwas}}/magma/magma_set_conditional.done"
  output:
    touch(f"{outdir}/reference/gwas_sumstat/{{gwas}}/magma/snplists/create_set_snplists.done")
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/create_set_snplists-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/create_set_snplists-{{gwas}}.log"
  params:
    config_file=config['config_file']
  shell:
    "Rscript ../Scripts/magma/set_extractor.R \
      --config {params.config_file} \
      --gwas {wildcards.gwas} > {log} 2>&1"
  
rule run_create_set_snplists:
  input: 
      lambda w: expand(f"{outdir}/reference/gwas_sumstat/{{gwas}}/magma/snplists/create_set_snplists.done", gwas=gwas_list_df_eur['name'])
  output: 
      touch(f"{outdir}/reference/gwas_sumstat/create_set_snplists_all_gwas.done")

# Create a file listing gwas with significant gene sets/set-specific SNP lists
checkpoint set_reporter:
  input:
    f"{outdir}/reference/gwas_sumstat/create_set_snplists_all_gwas.done"
  output:
    f"{outdir}/reference/gwas_sumstat/set_reporter.txt"
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/set_reporter.txt"
  log:
    f"{outdir}/reference/logs/set_reporter.log"
  params:
    config_file=config['config_file']
  shell:
    "Rscript ../Scripts/magma/set_reporter.R \
      --config {params.config_file}"
      
########
# Calculate stratified PGS
########

# Prepare score files for stratified PGS
rule pgs_stratifier:
  input:
    f"{outdir}/reference/gwas_sumstat/set_reporter.txt",
    rules.prep_pgs.input
  threads: config['cores_prep_pgs']
  output:
    touch(f"{outdir}/reference/pgs_score_files/pgs_stratifier.done"),
  conda:
    "../envs/analysis.yaml"
  benchmark:
    f"{outdir}/reference/benchmarks/pgs_stratifier.txt"
  log:
    f"{outdir}/reference/logs/pgs_stratifier.log"
  params:
    testing=config["testing"],
    config_file = config["config_file"]
  shell:
    "Rscript ../Scripts/pgs_methods/pgs_stratifier.R \
      --config {params.config_file} \
      --plink2 plink2 \
      --test {params.testing} \
      --n_cores {threads} > {log} 2>&1"

# Target sample scoring
rule target_pgs_partitioned_i:
  resources:
    mem_mb=config['mem_target_pgs'],
    time_min=1000
  threads: config['cores_target_pgs']
  input:
    f"{outdir}/reference/target_checks/{{name}}/ancestry_reporter.done",
    f"{outdir}/reference/gwas_sumstat/set_reporter.txt",
    lambda w: f"{outdir}/reference/target_checks/{{name}}/pc_projection-TRANS.done" if w.population == "TRANS" else [],
    f"{outdir}/reference/pgs_score_files/pgs_stratifier.done"
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/target_pgs_partitioned-{{population}}.done")
  benchmark:
    f"{outdir}/reference/benchmarks/target_pgs_partitioned_i-{{name}}-{{population}}.txt"
  log:
    f"{outdir}/reference/logs/target_pgs_partitioned_i-{{name}}-{{population}}.log"
  conda:
    "../envs/analysis.yaml"
  params:
    testing=config["testing"],
    config_file = config["config_file"]
  shell:
    "Rscript ../Scripts/target_scoring/target_scoring_partitioned_pipeline.R \
      --config {params.config_file} \
      --name {wildcards.name} \
      --population {wildcards.population} \
      --plink2 plink2 \
      --test {params.testing} \
      --n_cores {threads} > {log} 2>&1"

rule target_pgs_partitioned_all:
  input:
    lambda w: expand(f"{outdir}/reference/target_checks/{{name}}/target_pgs_partitioned-{{population}}.done", name=w.name, population = ancestry_munge(w.name, scaling = config["pgs_scaling"]))
  output:
    touch(f"{outdir}/reference/target_checks/{{name}}/target_pgs_partitioned.done")

rule target_pgs_partitioned:
  input:
    expand(f"{outdir}/reference/target_checks/{{name}}/target_pgs_partitioned.done", name=target_list_df['name'])
