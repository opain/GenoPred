# Focused pipeline checks

Run `Rscript pipeline/tests/scoring_regression.R` and
`python3 pipeline/tests/test_target_paths.py` from the repository root.
`test_workflow_contract.py` additionally requires Snakemake 7 on `PATH` (or set
`SNAKEMAKE` to its executable). The R integration scripts require the pipeline
analysis environment; they use temporary fixtures and mocked external tools.

`benchmark_score_io.R <score.gz> <ref.chr-prefix> <scratch-dir>` compares the
old per-chromosome extraction with one-pass sharding and checks that their
outputs are identical. It writes and removes temporary files in the supplied
scratch directory. On the mounted shared filesystem, one warm-cache sequential
run on PGS002804 (1,204,449 reference variants) took 14.9 s for the old
extraction including output comparisons, versus 5.0 s for sharding plus 1.9 s
for variant validation. The 10,319,088-byte compressed score produced
13,050,366 bytes of temporary weight shards; the reference map shared by all
scores added 19,993,438 bytes per scoring invocation. Timings and storage depend on score
column count, filesystem caching, concurrency and temporary-directory placement.
