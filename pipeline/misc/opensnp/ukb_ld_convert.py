#!/usr/bin/env python3
"""
ukb_ld_convert.py — Round 9 Section 2 converter.

For chr22 across all six UKB populations, read each population's Zarr store
from Zenodo 14614207 and emit one .npz per canonical block. R side then
converts .npz → .rds.

Reconciliation: the 6 populations use different per-pop LDetect maps
(EUR 24, EAS 20, AFR 34, CSA 20, AMR 34, MID 24 chr22 blocks; only 2–3
break points shared with EUR). To keep the block structure common across
populations while never fabricating cross-native-block correlations, we
take the **common refinement** of the 6 partitions — union of all
breakpoints. Every canonical sub-block is fully contained in each
population's native LDetect block, so all needed correlations are stored.

Storage layout (per pop): CSR upper-triangle, no diagonal.
Row i has (native_block_end - i - 1) entries, correlations to SNPs
i+1, i+2, ..., native_block_end-1. Reconstructed r = int8 / 127.

Alignment: canonical alleles = HM3 pvar (ref, alt). For each pop's row,
if UKB a1==alt and a2==ref (HM3 orientation) → no flip. If reversed → flip
both rows/cols and use 1-af1 as the a1-frequency. Any mismatch → drop.

Output: pipeline/misc/opensnp/meld_ld_ukb/chr22/_npz/block_XX.npz
"""
import sys, csv
from pathlib import Path
import numpy as np
import zarr

MISC = Path('/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp')
UKB  = Path('/users/k1806347/oliverpainfel/Data/ukb/zenodo_14614207')
OUT  = MISC / 'meld_ld_ukb' / 'chr22'
NPZ  = OUT / '_npz'
NPZ.mkdir(parents=True, exist_ok=True)

POPS = ['EUR','EAS','AFR','CSA','AMR','MID']

HM3_PVAR = Path('/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ref/ref.chr22.pvar')


def read_pop(pop):
    """Read one population's chr22 Zarr store."""
    root = UKB / pop / 'chr_22'
    z = zarr.open(str(root), mode='r')
    attrs = dict(z.attrs)
    m = z['metadata']
    md = m['maf'][:]
    d = dict(
        pop=pop,
        a1=np.asarray(m['a1'][:], dtype=object),
        a2=np.asarray(m['a2'][:], dtype=object),
        bp=np.asarray(m['bp'][:], dtype=int),
        cm=np.asarray(m['cm'][:], dtype=float),
        af1=np.asarray(md, dtype=float),
        snps=np.asarray(m['snps'][:], dtype=object),
        data=np.asarray(z['matrix/data'][:], dtype=np.int8),
        indptr=np.asarray(z['matrix/indptr'][:], dtype=np.int64),
        blocks=[tuple(b) for b in attrs['Estimator properties']['LD blocks']],
        build=attrs.get('Genome build', 'unknown'),
        n_samples=int(attrs.get('Sample size', -1)),
    )
    # Native block SNP-index bounds: for each native block, find the SNP-index
    # range [ns, ne) that falls within it. This is what indptr layout is keyed
    # to — indptr[i+1] - indptr[i] = ne - i - 1 for i in [ns, ne).
    d['native_snp_bounds'] = _native_snp_bounds(d['bp'], d['blocks'], d['indptr'])
    return d


def _native_snp_bounds(bp, blocks, indptr):
    """For each native block (bp_start, bp_end), return (snp_ns, snp_ne) —
    the SNP-index range within the flat variant list. Verify via indptr
    layout: for i in [ns, ne), indptr[i+1] - indptr[i] should equal ne-i-1."""
    out = []
    for (bs, be) in blocks:
        ns = int(np.searchsorted(bp, bs, side='left'))
        ne = int(np.searchsorted(bp, be, side='left'))
        # Verify layout
        if ne > ns:
            widths = indptr[ns+1:ne+1] - indptr[ns:ne]
            expected = np.arange(ne - ns - 1, -1, -1)
            if not np.array_equal(widths, expected):
                raise RuntimeError(
                    f'block bp=[{bs},{be}) snp=[{ns},{ne}): row widths do not '
                    f'match upper-tri layout. widths[:5]={widths[:5]}, '
                    f'expected[:5]={expected[:5]}'
                )
        out.append((ns, ne))
    return out


def get_r(pop_data, i, j):
    """Return r(i, j) for two SNPs in the same native block, aligned to
    the population's own a1 convention."""
    if i == j:
        return 1.0
    if i > j:
        i, j = j, i
    offset = j - i - 1
    return pop_data['data'][pop_data['indptr'][i] + offset] / 127.0


def build_canon_block_R(pop_data, snp_indices):
    """Given a canonical sub-block's SNP indices for this population (must all
    be within a single native block), reconstruct the M x M dense R matrix."""
    M = len(snp_indices)
    R = np.zeros((M, M), dtype=np.float32)
    for a in range(M):
        R[a, a] = 1.0
        i = snp_indices[a]
        # Row-i entries at indptr[i]:indptr[i+1] are correlations to SNPs
        # i+1, ..., native_block_end - 1.
        row = pop_data['data'][pop_data['indptr'][i]:pop_data['indptr'][i+1]]
        for b in range(a+1, M):
            j = snp_indices[b]
            r = row[j - i - 1] / 127.0
            R[a, b] = r
            R[b, a] = r
    return R


def common_refinement(all_pop_blocks):
    """Union of all breakpoints across populations → canonical sub-block
    boundaries."""
    breaks = set()
    for blocks in all_pop_blocks:
        for (s, e) in blocks:
            breaks.add(s)
            breaks.add(e)
    sb = sorted(breaks)
    return [(sb[i], sb[i+1]) for i in range(len(sb)-1)]


def load_hm3():
    """rsid → (ref, alt) from HM3 chr22 pvar."""
    triples = {}
    with open(HM3_PVAR) as f:
        for row in csv.reader(f, delimiter='\t'):
            if row and row[0].startswith('#'):
                continue
            chrom, pos, rsid, ref, alt = row[0], row[1], row[2], row[3], row[4]
            if len(ref) == 1 and len(alt) == 1:
                triples[rsid] = (ref, alt)
    return triples


def native_block_for_bp(native_bounds_bp, bp):
    """Return the index of the native block containing bp, else -1."""
    for k, (s, e) in enumerate(native_bounds_bp):
        if s <= bp < e:
            return k
    return -1


def main():
    hm3 = load_hm3()
    print(f'HM3 chr22 pvar (SNP-only): {len(hm3)} rsids')

    pop_data = {}
    print('\nLoading populations...')
    for p in POPS:
        d = read_pop(p)
        pop_data[p] = d
        print(f'  {p}: n_snps={len(d["snps"])}, N={d["n_samples"]}, '
              f'blocks={len(d["blocks"])}, build={d["build"]}')

    # Verify all builds are b37 and that block ranges cover the same span.
    for p in POPS:
        assert pop_data[p]['build'] == 'GRCh37', f'{p} build != GRCh37'
    print('  all populations on GRCh37 ✓')

    # Common refinement of block partitions.
    canon = common_refinement([pop_data[p]['blocks'] for p in POPS])
    print(f'\ncommon refinement: {len(canon)} sub-blocks')
    print(f'  first 3: {canon[:3]}')
    print(f'  last 3:  {canon[-3:]}')

    # Per-pop rsid → global-snp-index map
    rsid_to_idx = {p: {rs: i for i, rs in enumerate(pop_data[p]['snps'])} for p in POPS}

    # Build per-canonical-block dataset.
    per_block_counts = []
    total_variants = 0

    for bi, (bs, be) in enumerate(canon):
        block_id = bi + 1

        # Candidate rsids = union of pop rsids in this bp range, then ∩ HM3.
        cand = set()
        for p in POPS:
            bp = pop_data[p]['bp']
            i0 = int(np.searchsorted(bp, bs, side='left'))
            i1 = int(np.searchsorted(bp, be, side='left'))
            cand.update(pop_data[p]['snps'][i0:i1].tolist())
        cand = [rs for rs in cand if rs in hm3]

        # Retain only rsids present in every population and with alleles
        # matching HM3 (direct or swapped).
        keep = []
        for rs in cand:
            ok = True
            for p in POPS:
                if rs not in rsid_to_idx[p]:
                    ok = False; break
                i_p = rsid_to_idx[p][rs]
                a1u, a2u = pop_data[p]['a1'][i_p], pop_data[p]['a2'][i_p]
                ref_h, alt_h = hm3[rs]
                if not ((a1u == alt_h and a2u == ref_h) or (a1u == ref_h and a2u == alt_h)):
                    ok = False; break
            if ok:
                keep.append(rs)

        if not keep:
            per_block_counts.append(dict(block_id=block_id, m=0, bp=[bs,be]))
            continue

        # Order by bp (use EUR bp — populations share position for shared rsid)
        eur_bp = np.array([pop_data['EUR']['bp'][rsid_to_idx['EUR'][rs]] for rs in keep])
        order = np.argsort(eur_bp)
        keep = [keep[o] for o in order]

        M = len(keep)
        total_variants += M

        # For each population, ensure all `keep` SNPs are in one native block
        # (by construction of common refinement, this must be true).
        pop_local_ix = {}
        for p in POPS:
            ixs = np.array([rsid_to_idx[p][rs] for rs in keep])
            pop_local_ix[p] = ixs
            # Verify all within a single native block
            native = [native_block_for_bp(
                [(s, e) for (s, e) in pop_data[p]['blocks']],
                int(pop_data[p]['bp'][i])) for i in ixs]
            if len(set(native)) != 1:
                raise RuntimeError(
                    f'block {block_id} pop {p}: SNPs span multiple native '
                    f'blocks: {set(native)}. Common-refinement invariant '
                    f'violated.')

        # Build per-pop R by iterating pairs (use the stored upper triangle).
        block_out = dict(
            chr=np.int32(22),
            block_id=np.int32(block_id),
            block_bp=np.array([bs, be], dtype=np.int64),
            SNP=np.array(keep, dtype=object),
            BP=np.array([pop_data['EUR']['bp'][rsid_to_idx['EUR'][rs]] for rs in keep], dtype=np.int32),
            cM=np.array([pop_data['EUR']['cm'][rsid_to_idx['EUR'][rs]] for rs in keep], dtype=np.float32),
            A1_canon=np.array([hm3[rs][1] for rs in keep], dtype=object),  # HM3 ALT
            A2_canon=np.array([hm3[rs][0] for rs in keep], dtype=object),  # HM3 REF
        )

        for p in POPS:
            ixs = pop_local_ix[p]
            flip = np.zeros(M, dtype=bool)
            f_a1_local = pop_data[p]['af1'][ixs].copy()
            for k, rs in enumerate(keep):
                a1u = pop_data[p]['a1'][ixs[k]]
                ref_h, alt_h = hm3[rs]
                if a1u == ref_h:  # UKB a1 == HM3 ref → flip
                    flip[k] = True
            # Build R for this pop using pairwise get_r on native storage.
            R = build_canon_block_R(pop_data[p], ixs)
            # Apply flip: R_ij → s_i s_j R_ij; f → 1-f for flipped SNPs.
            if flip.any():
                D = np.where(flip, -1.0, 1.0).astype(np.float32)
                R = R * np.outer(D, D)
                f_a1_local = np.where(flip, 1.0 - f_a1_local, f_a1_local)
            np.fill_diagonal(R, 1.0)
            v = 2.0 * f_a1_local * (1.0 - f_a1_local)
            block_out[f'R_{p}'] = R.astype(np.float32)
            block_out[f'f_{p}'] = f_a1_local.astype(np.float32)
            block_out[f'v_{p}'] = v.astype(np.float32)

        np.savez(NPZ / f'block_{block_id:03d}.npz', **block_out)
        per_block_counts.append(dict(block_id=block_id, m=M, bp=[bs, be]))
        if block_id % 10 == 0:
            print(f'  wrote block {block_id}/{len(canon)}  M={M}')

    # Summary
    print(f'\nsummary:')
    print(f'  canonical sub-blocks: {len(canon)}')
    print(f'  non-empty sub-blocks: {sum(1 for r in per_block_counts if r["m"] > 0)}')
    print(f'  total variants across sub-blocks: {total_variants}')
    small = [r for r in per_block_counts if 0 < r["m"] < 30]
    print(f'  sub-blocks with 1..29 variants (unusable in M2): {len(small)}')
    print(f'  sub-blocks with 30+ variants (M2-usable): '
          f'{sum(1 for r in per_block_counts if r["m"] >= 30)}')
    with open(OUT / 'block_counts.csv', 'w') as fh:
        fh.write('block_id,m,bp_start,bp_end\n')
        for r in per_block_counts:
            fh.write(f'{r["block_id"]},{r["m"]},{r["bp"][0]},{r["bp"][1]}\n')
    print(f'wrote {OUT / "block_counts.csv"}')


if __name__ == '__main__':
    main()
