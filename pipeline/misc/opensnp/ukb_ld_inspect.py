#!/usr/bin/env python3
"""
ukb_ld_inspect.py — Round 9 Section 2 inspection.

Report the on-disk structure of the extracted UKB LD matrices from Zenodo
14614207 (Zabad et al.). This runs BEFORE writing the converter and its
output is what gets brought back to the user for confirmation.

Reports:
  - directory tree under one population, one chromosome
  - variant-table schema and a few rows
  - how blocks are keyed
  - allele-frequency column and its range
  - genome build
  - int8 -> correlation conversion (scale, whether diag is stored, min/max)
  - LDetect block-boundary match with Round 8's fourier_ls-chr22.bed
"""

import os, sys, json
from pathlib import Path
import numpy as np
import zarr

BASE = Path('/users/k1806347/oliverpainfel/Data/ukb/zenodo_14614207')
POPS = ['EUR','CSA','AFR','EAS','MID','AMR']
R8_BED = Path('/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ld_blocks/EUR/fourier_ls-chr22.bed')


def hdr(t):
    print(f'\n=== {t} ===')

# ---------------------------------------------------------------------------
hdr('1. Directory tree — EUR chr_22')
eur22 = BASE / 'EUR' / 'chr_22'
for p in sorted(eur22.rglob('*'))[:60]:
    rel = p.relative_to(eur22)
    depth = len(rel.parts)
    print(f"  {'  '*(depth-1)}{rel.parts[-1]}  {'(dir)' if p.is_dir() else f'{p.stat().st_size}b'}")

# ---------------------------------------------------------------------------
hdr('2. Zarr group inspection — EUR chr_22')
# Open the top-level Zarr group.
try:
    z = zarr.open(str(eur22), mode='r')
    print(f'root type: {type(z).__name__}')
    print(f'root attrs: {dict(z.attrs)}')
    print('tree:')
    print(z.tree())
except Exception as e:
    print(f'zarr.open failed: {e}')
    # fall back to reading each array we can find
    for f in sorted(eur22.rglob('.zarray')):
        with open(f) as h:
            spec = json.load(h)
        print(f'  {f.parent.relative_to(eur22)}: shape={spec.get("shape")} '
              f'dtype={spec.get("dtype")} chunks={spec.get("chunks")}')

# ---------------------------------------------------------------------------
hdr('3. Variant-table schema')
# The metadata/ subdir typically holds per-variant arrays keyed by field.
meta = eur22 / 'metadata'
if meta.is_dir():
    for sub in sorted(meta.iterdir()):
        if sub.is_dir():
            zar = sub / '.zarray'
            if zar.exists():
                with open(zar) as h:
                    spec = json.load(h)
                arr = zarr.open(str(sub), mode='r')
                sample = np.asarray(arr[:5])
                print(f"  {sub.name}: dtype={spec['dtype']}  shape={spec['shape']}  sample={sample}")

# ---------------------------------------------------------------------------
hdr('4. Block layout')
# The matrix/ subdir. Blocks may be keyed by index or by (start,end) — look.
mat = eur22 / 'matrix'
if mat.is_dir():
    print('matrix/ contents:')
    for sub in sorted(mat.iterdir())[:8]:
        print(f'  {sub.name}  {"(dir)" if sub.is_dir() else sub.stat().st_size}')
    # Report .zattrs at matrix root and one block
    if (mat / '.zattrs').exists():
        with open(mat / '.zattrs') as h:
            print(f"matrix .zattrs: {h.read()}")
    # First array
    zar_candidates = list(mat.rglob('.zarray'))
    if zar_candidates:
        for z1 in zar_candidates[:3]:
            with open(z1) as h:
                spec = json.load(h)
            print(f'  {z1.parent.relative_to(mat)}: dtype={spec["dtype"]}  shape={spec["shape"]}  '
                  f'chunks={spec.get("chunks")}  compressor={spec.get("compressor")}')

# ---------------------------------------------------------------------------
hdr('5. magenpy-native load')
# Try to load via magenpy to get the block structure and int8 -> r scaling.
try:
    import magenpy as mp
    print(f'magenpy version: {mp.__version__}')
    # LDMatrix loader — the Zenodo docs indicate one Zarr store per chromosome
    ld = mp.LDMatrix.from_path(str(eur22))
    print(f'type: {type(ld).__name__}')
    print(f'n_snps: {getattr(ld, "n_snps", "?")}')
    print(f'estimator: {getattr(ld, "ld_estimator", "?")}')
    print(f'has store_attrs: {list(ld.store_attrs) if hasattr(ld, "store_attrs") else "?"}')
    if hasattr(ld, 'attrs'):
        print(f'attrs: {dict(ld.attrs)}')
    # Block-keying diagnostic
    for attr in ('ld_boundaries', 'ld_block_bounds', 'blocks', 'window_size'):
        if hasattr(ld, attr):
            v = getattr(ld, attr)
            print(f'  {attr}: type={type(v).__name__}, shape/len={getattr(v,"shape",len(v) if hasattr(v,"__len__") else "?")}')
            if attr == 'ld_boundaries':
                print(f'  first 3 rows: {np.asarray(v)[:3]}')
    # Sample the first block's correlations
    try:
        r0 = ld.load_data()  # or similar; interface may vary by version
        print(f'load_data returned: type={type(r0).__name__}')
    except Exception as ee:
        print(f'load_data not directly usable: {ee}')
except ImportError:
    print('magenpy not importable')

# ---------------------------------------------------------------------------
hdr('6. int8 -> r conversion inspection')
# Read one raw int8 chunk to see values, without magenpy.
# Look for the data array under matrix/
data_dir = mat / 'data' if (mat / 'data').is_dir() else mat
if data_dir.is_dir():
    zar = data_dir / '.zarray'
    if zar.exists():
        with open(zar) as h:
            spec = json.load(h)
        print(f'data .zarray: {spec}')
        # Load first chunk
        arr = zarr.open(str(data_dir), mode='r')
        raw = np.asarray(arr[:100])  # first 100 entries
        print(f'first 100 raw values: min={raw.min()}, max={raw.max()}, dtype={raw.dtype}')
        # The typical convention is int8 with scale 1/127 for correlations in [-1,1].
        # If .zattrs on the array specifies a scale factor, use it; otherwise this is
        # a working hypothesis to be confirmed via magenpy.
        for attrs_path in (data_dir / '.zattrs', mat / '.zattrs', eur22 / '.zattrs'):
            if attrs_path.exists():
                print(f'\n{attrs_path.relative_to(eur22)} contents:')
                print(attrs_path.read_text())

# ---------------------------------------------------------------------------
hdr('7. LDetect block boundaries comparison')
# Compare UKB block starts/ends (from ld_boundaries) against Round 8's BED file.
r8_bounds = np.loadtxt(R8_BED, skiprows=1, dtype=int, usecols=(1,2))
print(f'Round 8 bed boundaries (chr22): {len(r8_bounds)} pairs')
print(f'  first 3: {r8_bounds[:3]}')
print(f'  last 3:  {r8_bounds[-3:]}')

try:
    import magenpy as mp
    ld = mp.LDMatrix.from_path(str(eur22))
    if hasattr(ld, 'ld_boundaries'):
        ukb_bounds = np.asarray(ld.ld_boundaries)
        print(f'UKB boundaries (chr22): shape={ukb_bounds.shape}')
        print(f'  first 3: {ukb_bounds[:3]}')
        print(f'  last 3:  {ukb_bounds[-3:]}')
except Exception as e:
    print(f'boundary comparison via magenpy failed: {e}')

# ---------------------------------------------------------------------------
hdr('8. AF availability and range')
# Look for an AF/MAF-like array.
for cand in ('maf', 'MAF', 'af', 'AF', 'allele_freq', 'freq'):
    p = meta / cand
    if p.is_dir() and (p / '.zarray').exists():
        arr = zarr.open(str(p), mode='r')
        v = np.asarray(arr[:])
        print(f'  {cand}: n={len(v)} min={v.min():.4f} max={v.max():.4f} '
              f'mean={v.mean():.4f} sample={v[:5]}')

print('\n=== DONE ===')
