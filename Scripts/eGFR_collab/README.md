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
│   ├── pseudovalidated.score.gz    PRIMARY: 1 PGS per method per source, no tuning needed
│   └── full_grid.score.gz          OPTIONAL: every hyperparameter (large; only for CV)
└── scripts/
    ├── 01_compute_pgs.sh           PLINK2 --score template
    └── 02_evaluate_pgs.R           per-PGS incremental R²; optional CV across methods
```

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

Edit `scripts/01_compute_pgs.sh` to point at your PLINK2 binary and target
`.pgen/.pvar/.psam` files, then run:

```bash
bash scripts/01_compute_pgs.sh
```

Output is `pgs_out/pseudovalidated.sscore` with one column per PGS column
(`<column_name>_AVG` is the per-individual average effect; `<column_name>_SUM`
is the sum).

> **Variant overlap.** PLINK2 prints how many variants matched. If far fewer
> than expected, double-check chromosome/build (the bundle is GRCh38) and SNP
> ID format.

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

### Default mode (recommended)

```bash
Rscript scripts/02_evaluate_pgs.R \
  --pgs        pgs_out/pseudovalidated.sscore \
  --pheno      egfr.tsv \
  --covar      covariates.tsv \
  --catalogue  catalogue.tsv \
  --target_pop EUR \
  --output     results_EUR
```

For each pseudo-validated PGS relevant to the target population, fits

```
phenotype ~ age + sex + PCs + scale(PGS)
```

and reports the **incremental R²** vs the covariates-only model with a
bootstrap 95% CI. No tuning, no model selection — the incremental R² is an
honest estimate of out-of-sample predictive utility.

Output: `<output>.per_pgs.tsv`.

### Optional CV mode

```bash
Rscript scripts/02_evaluate_pgs.R ... --cv --folds 5
```

Runs nested k-fold CV across **all** PGS columns to estimate an unbiased
"best PGS across methods" incremental R². Inner folds pick the best column on
training data; outer folds score it on held-out data. Without nested CV,
picking the highest R² across hundreds of columns dramatically inflates the
estimate — a common pitfall when moving from PRSice-style sweeping to many
methods at once.

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
