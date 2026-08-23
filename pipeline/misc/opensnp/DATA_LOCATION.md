# MELD data, results, and LD panels — where they live

The scripts in this directory analyse ~30 GB of LD reference matrices,
~10 MB of harmonised sumstats, ~10 MB of result CSVs, and rendered
notebook HTMLs. Those live **outside the git repo** at:

```
/users/k1806347/oliverpainfel/Analyses/MELD/
├── README.md
├── data/                         inputs derived from external sources
│   ├── gbmi_r8_harmonised/       7 GBMI-trait chr22 RDSes + per-pop AF cache
│   ├── r7bS4/close_snp_sets.rds  precomputed intersections for R7bS4 M5
│   └── manifest_GBMI_summary_statistics.csv
├── ld/                           LD reference matrices (~30 GB)
│   ├── 1kg_hgdp/                 24 chr22 block RDSes × 5 pops (was meld_ld/)
│   ├── 1kg_hgdp_on_r9blocks/     R8 LD on the R9 sub-block structure
│   ├── ukb/                      67 sub-blocks converted from Zenodo 14614207
│   ├── ref_empirical/            1KG+HGDP empirical AF/var per pop
│   └── private_ldpred2_ldref/    AMR-slot symlink hijack wrappers *
├── results/                      one subdirectory per round
│   ├── r4/  r5/  r6/  r6b/  r7/  r7b/  r8/  r9/  r10/
└── rendered/                     8 notebook HTMLs (r6, r6b, r7 §1-2, r7bA,
                                  r7bS4, r8, r9, r10)
```

Every MELD script sources `meld_paths.R` (in this directory) to pick up
the canonical paths, so relocating the analysis tree only requires
editing `MELD_ROOT` in that one file.

## Common questions

**Where is trait X's harmonised chr22 sumstats?**
`~/Analyses/MELD/data/gbmi_r8_harmonised/gbmi_r8_<TRAIT>_chr22.rds`.
7 traits: Asthma, COPD, Gout, HF, IPF, Stroke, VTE.

**Where is Round N's M2 output?**
`~/Analyses/MELD/results/rN/` — the CSV names are unchanged from the
original notebooks; only the directory changed.

**Where are the block LD RDSes for chr22?**
`~/Analyses/MELD/ld/1kg_hgdp/chr22/block_XXXX.rds`. UKB variants at
`~/Analyses/MELD/ld/ukb/chr22/`.

**Where are the rendered notebooks?**
`~/Analyses/MELD/rendered/opensnp_meld_rN.html`. The source `.Rmd`
stays in `docs/` in this repo.

## \* about `private_ldpred2_ldref_*`

The seven `private_ldpred2_ldref_<label>/AMR/` symlinks were originally
placed under `pipeline/misc/opensnp/` because the GenoPred pipeline
rules look for them there (they hijack the AMR slot for the Round 7b
Section 4 LDpred2 refdir experiment). After the relocation, these
symlinks live under `~/Analyses/MELD/ld/private_ldpred2_ldref/`. **If
the pipeline is re-invoked from this checkout to run any R7bS4
LDpred2 config, the seven `private_ldpred2_ldref_<label>` symlinks
must be re-created under `pipeline/misc/opensnp/`** — the pipeline
does not look up the new location. `misc/opensnp/setup_r7bS4_ldpred2.sh`
still writes them into the pipeline directory; re-running it after
setting `BUILT=~/Analyses/MELD/ld/ldpred2_ref_meld` restores the setup.

## To relocate the whole analysis tree

Edit `MELD_ROOT` in `meld_paths.R` (single line change). All scripts
pick up the new location on their next run.
