# eGFR PGS share bundle

This bundle contains polygenic score (PGS) weights derived from two eGFR GWAS
(EUR and AFR) using multiple PGS methods, plus a small set of scripts for
calculating and evaluating the PGS in your target sample (e.g. *All of Us*).

The PGS weights were generated with [GenoPred](https://opain.github.io/GenoPred/)
and have already been QC'd: ambiguous variants removed, variants aligned to the
positive strand, weights expressed for the A1 allele. Almost all variants
should match your target genotype data — please report the overlap PLINK2
prints (`--score: N variants processed`) as a sanity check.

## Contents

```
egfr_pgs_share/
├── README.md                       (this file)
├── catalogue.tsv                   (one row per PGS column: method, source, tuning info)
├── scores/
│   ├── ptclump_only.score.gz       SENSITIVITY: p+T at every P-threshold, PRSice-style
│   ├── pseudovalidated.score.gz    SumStatTune: 1 pseudo-val PGS per method per source
│   └── full_grid.score.gz          IndivTune: every hyperparameter (large; for CV tuning)
└── scripts/
    ├── 01_compute_pgs.sh           PLINK2 --score template
    ├── 02_evaluate_pgs.R           SumStatTune + IndivTune + ptclump sensitivity
    └── 03_plot.R                   figure: R per (method × source × target)
```

## Three analyses

1. **ptclump sensitivity** — every P-threshold for pT+clump, for direct
   comparison against your existing PRSice results.
2. **SumStatTune** — for each (method × source GWAS × target ancestry), the R
   of the **pseudo-validated** PGS (no target tuning, derived from summary
   statistics alone). LEOPARD-combined columns provide the multi-source
   (EUR+AFR) row.
3. **IndivTune** — for each (method × source GWAS × target ancestry), the R
   from **nested k-fold CV in the target sample**. Single-source: inner CV
   picks the best hyperparameter. EUR+AFR: independent inner CV picks the
   best per-source column, outer fold joint-fits both.

Run each per target ancestry (with `--keep`), then plot SumStat and Indiv R
side-by-side with `03_plot.R`.

## Score file format

Plain text, space-separated, gzipped. Columns:

```
SNP  A1  A2  <pgs_col_1>  <pgs_col_2>  ...
```

Each `<pgs_col_*>` is a set of SNP weights (β for the A1 allele) for one
method / source GWAS / hyperparameter combination. Column names are unique
and match `catalogue.tsv$column_name`.

## The catalogue

`catalogue.tsv` describes every PGS column. Key fields:

| Field                | Meaning                                                                   |
| -------------------- | ------------------------------------------------------------------------- |
| `column_name`        | Header of the corresponding column in the score file                      |
| `method`             | PGS method (`ldpred2`, `prscs`, `prscsx`, `quickprs_multi`, etc.)         |
| `source_gwas`        | Source GWAS used to derive the weights (`egfr_eur`, `egfr_afr`, `egfr`)   |
| `source_population`  | Population of the source GWAS (`EUR`, `AFR`, `EUR+AFR`, `META`)           |
| `hyperparameter`     | Method-specific parameter (e.g. `beta_auto`, `phi_1e-04`)                 |
| `is_pseudovalidated` | `TRUE` if this column requires **no tuning** — recommended for headline use |
| `target_population`  | `EUR`/`AFR` if the column was tuned **for a specific target population**; `NA` otherwise |

The pseudo-validated subset spans:

* `ldpred2` (`beta_auto`), `prscs` (`phi_auto`), `sbayesrc`, `quickprs`,
  `dbslmm` — single-source, auto-tuned. Use the one matching your target
  population, e.g. `egfr_eur` source for an EUR target.
* `prscsx` (`META_phi_auto`) — multi-source, cross-population pseudo-val.
* `quickprs_multi`, `prscs_multi`, `megaprs_multi`, `ptclump_multi` — LEOPARD
  weighted combinations of population-specific PGS, **tuned for a specific
  target population** (see `target_population` column).

Methods without a pseudo-validated single model (e.g. `ptclump`, `lassosum`,
`megaprs` hyperparameter grids) are only included in `full_grid.score.gz` and
should be used via CV-tuning (see below).

## Step 1 — Compute PGS with PLINK2

Edit `scripts/01_compute_pgs.sh` and set **one** of:

* `TARGET_PFILE=/path/to/all_chr_data` if you have a single whole-genome
  PLINK2 `.pgen/.pvar/.psam` set, **or**
* `TARGET_PFILE_CHR=/path/to/chr{CHR}/genotypes` if your data is split per
  chromosome (the `{CHR}` placeholder is substituted with 1, 2, …, 22). The
  script will score each chromosome separately and sum the contributions to
  produce a single `.sscore`.

Then run:

```bash
bash scripts/01_compute_pgs.sh
```

Output is `pgs_out/pseudovalidated.sscore` with one column per PGS column
(`<column_name>_AVG` is the per-individual average effect over scored
variants; `<column_name>_SUM` is the dosage-weighted sum).

> **Variant overlap.** PLINK2 prints how many variants matched per
> chromosome. Summed across chromosomes should be close to the total in the
> score file (~1.2M for this bundle). If far fewer, double-check
> chromosome/build (the bundle is GRCh38) and SNP ID format.

## Step 2 — Evaluate PGS

`scripts/02_evaluate_pgs.R` is for continuous outcomes only (eGFR or a
transformation of your choosing — apply the transformation before passing the
phenotype in).

### Inputs

* **`--pgs`** PLINK2 `.sscore` from Step 1.
* **`--pheno`** TSV with columns `IID <pheno>` (or `FID IID <pheno>`).
* **`--covar`** TSV with `IID age sex PC1 PC2 ...` (or with `FID`). Sex
  numeric, PCs as you normally use them (10–20 typical).
* **`--catalogue`** `catalogue.tsv` from this bundle.
* **`--target_pop`** `EUR` or `AFR` — filters which PGS columns are relevant.
* **`--keep`** (**strongly recommended**) PLINK keep file (`FID IID`,
  whitespace-separated, header optional) restricting the analysis to a single
  target population. PGS evaluation should be ancestry-stratified — without
  this flag the script runs on the full merged sample and emits a warning.

### Step 2a — ptclump sensitivity vs PRSice (recommended first)

If you've previously analysed these summary statistics with PRSice, run this
first to confirm GenoPred weights give comparable signal. Differences are
expected because GenoPred restricts to HapMap3 variants whereas PRSice uses
all GWAS×target overlap, but rank ordering of P-thresholds should be similar.

```bash
Rscript scripts/02_evaluate_pgs.R \
  --pgs        pgs_out/ptclump_only.sscore \
  --pheno      egfr.tsv \
  --covar      covariates.tsv \
  --catalogue  catalogue.tsv \
  --target_pop EUR \
  --keep       samples_EUR.keep \
  --output     results_EUR_ptclump \
  --all_columns
```

Output: `<output>.all_columns.tsv` with one row per P-threshold per source GWAS.
The maximum R picked across thresholds on a single sample is overfit — use
this only for comparison-to-PRSice purposes, not as a headline R.

### Step 2b — SumStatTune (headline)

```bash
Rscript scripts/02_evaluate_pgs.R \
  --pgs        pgs_out/pseudovalidated.sscore \
  --pheno      egfr.tsv \
  --covar      covariates.tsv \
  --catalogue  catalogue.tsv \
  --target_pop EUR \
  --keep       samples_EUR.keep \
  --output     results_EUR
```

Output: `<output>.sumstat_tune.tsv` with one row per (method, source, target),
giving R and SE for the pseudo-validated PGS. No target tuning — coefficients
are derived entirely from the summary statistics. Repeat for `--target_pop AFR
--keep samples_AFR.keep --output results_AFR`.

### Step 2c — IndivTune (CV in target sample)

```bash
Rscript scripts/02_evaluate_pgs.R \
  --pgs        pgs_out/full_grid.sscore \
  --pheno      egfr.tsv \
  --covar      covariates.tsv \
  --catalogue  catalogue.tsv \
  --target_pop EUR \
  --keep       samples_EUR.keep \
  --output     results_EUR \
  --indiv_tune --folds 5
```

Output: `<output>.indiv_tune.tsv` with the same shape, but R is from nested
k-fold CV (inner CV picks best hyperparameter / per-source column on training,
outer fold scores on held-out data). Joint EUR+AFR rows fit both sources
together in the outer fold.

### Step 2d — Plot

```bash
Rscript scripts/03_plot.R \
  --inputs "EUR=results_EUR,AFR=results_AFR" \
  --output figure_eGFR.png \
  --title  "eGFR PGS evaluation"
```

Renders a facet grid (rows = target, cols = source GWAS) with both
SumStatTune and IndivTune points + SE bars per method.

## Multi-ancestry analysis advice

Run the analysis **stratified by target ancestry** (separate EUR and AFR
samples). For each ancestry stratum:

1. **Single-source pseudo-validated PGS** (`is_pseudovalidated=TRUE`,
   `target_population=NA`): try the source matching the target ancestry first
   (e.g. `source_population=EUR` for an EUR target). Also try the cross-ancestry
   source — it is sometimes competitive.
2. **Cross-population pseudo-validated PGS**: include `prscsx` `META_phi_auto`
   (`source_population=META`). This combines EUR+AFR information at
   weight-fitting time and works regardless of target.
3. **Target-tuned multi-source PGS** (`target_population=<your target>`): the
   LEOPARD-combined `*_multi` columns are weighted combinations of
   population-specific PGS optimised for a specific target population. Use the
   column whose `target_population` matches your target stratum.

A weighted combination of EUR- and AFR-source PGS, with weights optimised for
the target population, is usually the best approach when you can't tune in the
target sample. The LEOPARD-combined columns deliver this without requiring
target-sample tuning. See [Pain et al. 2025] for details:
[GenoPred multi-ancestry PGS paper](https://opain.github.io/GenoPred/).

## Questions

Contact: oliver.pain@kcl.ac.uk
