# Automatic PUMAS-EN integration for one EUR novel GWAS. Synthetic fold
# summary statistics are exposed through gwas_list_df by dependencies.smk, so
# each selected GenoPred method is trained using its existing individual rule.

# Setup target materializes the dedicated rule environment on a compute node
# before PUMAS-EN is enabled for web submissions.
rule prepare_pumas_en_environment:
  input:
    env=f"{workflow.basedir}/envs/pumas.yaml",
    cpp=f"{pumas_en_code}/evaluation/CoordDescent.cpp"
  output:
    touch(f"{resdir}/software/PUMAS/pumas-environment.ready")
  conda:
    "../envs/pumas.yaml"
  shell:
    """
    Rscript -e 'stopifnot(all(vapply(c("data.table", "R.utils", "optparse", "BEDMatrix", "Rcpp", "jsonlite"), requireNamespace, logical(1), quietly=TRUE))); Rcpp::sourceCpp("{input.cpp}")' &&
      mkdir -p {resdir}/software/PUMAS && touch {output}
    """

pumas_en_train_outputs = [
  f"{outdir}/reference/pumas_en/{{gwas}}/upstream/{{gwas}}.gwas.omnibus.ite{fold}.txt"
  for fold in range(1, pumas_en_folds + 1)
]
pumas_en_xty_outputs = [
  f"{outdir}/reference/pumas_en/{{gwas}}/upstream/{{gwas}}.xty.omnibus.ite{fold}.txt"
  for fold in range(1, pumas_en_folds + 1)
]

rule pumas_en_subsample_i:
  resources:
    mem_mb=32000,
    time_min=1440
  threads: config['pumas_en_cores']
  input:
    cleaned=f"{outdir}/reference/gwas_sumstat/{{gwas}}/{{gwas}}-cleaned.gz",
    ld_blocks=pumas_en_ld_blocks,
    puma_subsample=f"{pumas_en_code}/PUMA-ensemble.subsampling.R",
    puma_commit=os.path.join(os.path.dirname(pumas_en_code), "COMMIT"),
    puma_manifest=pumas_en_manifest,
    sub_helpers=f"{pumas_en_code}/subsampling/helpers.R",
    eval_helpers=f"{pumas_en_code}/evaluation/helpers.R",
    adapter=f"{workflow.basedir}/../Scripts/pgs_methods/pumas_en_subsample.R",
    helpers=f"{workflow.basedir}/../Scripts/pgs_methods/pumas_en_helpers.R"
  output:
    folds=pumas_en_train_outputs,
    xty=pumas_en_xty_outputs,
    stats=f"{outdir}/reference/pumas_en/{{gwas}}/upstream/{{gwas}}.omnibus.forEVAL.txt"
  benchmark:
    f"{outdir}/reference/benchmarks/pumas_en_subsample_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/pumas_en_subsample_i-{{gwas}}.log"
  conda:
    "../envs/pumas.yaml"
  params:
    work_dir=lambda w: f"{outdir}/reference/pumas_en/{w.gwas}",
    code_dir=pumas_en_code,
    ld_blocks=pumas_en_ld_blocks,
    folds=pumas_en_folds,
    partitions=",".join(map(str, pumas_en_partitions))
  shell:
    "Rscript {input.adapter} "
    "--sumstats {input.cleaned} "
    "--gwas {wildcards.gwas} "
    "--work-dir {params.work_dir} "
    "--code-dir {params.code_dir} "
    "--ld-blocks {params.ld_blocks} "
    "--helpers {input.helpers} "
    "--folds {params.folds} "
    "--partitions {params.partitions} "
    "--threads {threads} > {log} 2>&1"


def pumas_en_full_score_paths(wildcards):
  return [
    f"{outdir}/reference/pgs_score_files/{method}/{wildcards.gwas}/ref-{wildcards.gwas}.score.gz"
    for method in pumas_en_methods
  ]


def pumas_en_fold_score_paths(wildcards):
  return [
    f"{outdir}/reference/pgs_score_files/{method}/{fold}/ref-{fold}.score.gz"
    for fold in pumas_fold_names
    for method in pumas_en_methods
  ]


rule prep_pgs_pumas_en_i:
  resources:
    mem_mb=64000,
    time_min=1440
  threads: config['pumas_en_cores']
  input:
    full_scores=pumas_en_full_score_paths,
    fold_scores=pumas_en_fold_score_paths,
    xty=lambda w: [
      f"{outdir}/reference/pumas_en/{w.gwas}/upstream/{w.gwas}.xty.omnibus.ite{fold}.txt"
      for fold in range(1, pumas_en_folds + 1)
    ],
    stats=lambda w: f"{outdir}/reference/pumas_en/{w.gwas}/upstream/{w.gwas}.omnibus.forEVAL.txt",
    ref_bed=f"{pumas_en_ref_prefix}.bed",
    ref_bim=f"{pumas_en_ref_prefix}.bim",
    ref_fam=f"{pumas_en_ref_prefix}.fam",
    puma_eval=f"{pumas_en_code}/PUMAS-ensemble.evaluation.R",
    puma_commit=os.path.join(os.path.dirname(pumas_en_code), "COMMIT"),
    puma_manifest=pumas_en_manifest,
    puma_model=f"{pumas_en_code}/evaluation/PUMAS-ensemble.R",
    puma_helpers=f"{pumas_en_code}/evaluation/helpers.R",
    puma_cpp=f"{pumas_en_code}/evaluation/CoordDescent.cpp",
    adapter=f"{workflow.basedir}/../Scripts/pgs_methods/pumas_en.R",
    helpers=f"{workflow.basedir}/../Scripts/pgs_methods/pumas_en_helpers.R"
  output:
    score=f"{outdir}/reference/pgs_score_files/pumas_en/{{gwas}}/ref-{{gwas}}.score.gz",
    provenance=f"{outdir}/reference/pgs_score_files/pumas_en/{{gwas}}/pumas_en.provenance.json"
  benchmark:
    f"{outdir}/reference/benchmarks/prep_pgs_pumas_en_i-{{gwas}}.txt"
  log:
    f"{outdir}/reference/logs/prep_pgs_pumas_en_i-{{gwas}}.log"
  conda:
    "../envs/pumas.yaml"
  params:
    methods=",".join(pumas_en_methods),
    full_scores=lambda w: ",".join(pumas_en_full_score_paths(w)),
    fold_scores=lambda w: ",".join(pumas_en_fold_score_paths(w)),
    work_dir=lambda w: f"{outdir}/reference/pumas_en/{w.gwas}",
    code_dir=pumas_en_code,
    ref_prefix=pumas_en_ref_prefix,
    folds=pumas_en_folds,
    partitions=",".join(map(str, pumas_en_partitions)),
    commit=pumas_en_commit,
    resource_manifest=pumas_en_manifest,
    testing=config["testing"]
  shell:
    "Rscript {input.adapter} "
    "--gwas {wildcards.gwas} "
    "--methods {params.methods} "
    "--full-scores {params.full_scores} "
    "--fold-scores {params.fold_scores} "
    "--folds {params.folds} "
    "--partitions {params.partitions} "
    "--threads {threads} "
    "--code-dir {params.code_dir} "
    "--work-dir {params.work_dir} "
    "--reference-prefix {params.ref_prefix} "
    "--ref-plink-chr {refdir}/ref.chr "
    "--test {params.testing} "
    "--out-score {output.score} "
    "--helpers {input.helpers} "
    "--pumas-commit {params.commit} "
    "--resource-manifest {params.resource_manifest} > {log} 2>&1"
