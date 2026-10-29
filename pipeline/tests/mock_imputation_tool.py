#!/usr/bin/env python3
"""Small deterministic command stub for the imputer integration test."""
import gzip
from pathlib import Path
import sys

args = sys.argv[1:]

def value(flag):
    return args[args.index(flag) + 1]

def write(path, content):
    Path(path).write_text(content)

ped = 'F I 0 0 1 -9 A A\n'
mapdata = '22 rs123 0 100\n'
if '--23file' in args:
    out = value('--out')
    write(out + '.ped', ped)
    write(out + '.map', mapdata)
elif '--file' in args:
    out = value('--out')
    write(out + '.ped', ped)
    write(out + '.map', mapdata)
elif '-check' in args:
    out = value('--output-log')
    write(out + '.snp.strand', 'type1 type2 pos main_id main_A main_B main_flippable ref_id ref_A ref_B ref_flippable\nStrand x 100 rs123 A A x rs123 A C x\n')
elif '--input-ped' in args:
    out = value('-O')
    write(out + '.haps', '22 rs123 100 A C 0 1 0 1\n')
    write(out + '.sample', 'ID_1 ID_2 missing\n0 0 0\nF I 0\n0 dummy 0\n')
elif '--hapsample2vcf' in args or '+fill-tags' in args:
    write(value('-o'), 'bcf\n')
elif 'index' in args:
    pass
elif '--h' in args:
    write(value('--o'), 'bcf\n')
elif '--gensample' in args:
    out = value('--gensample')
    with gzip.open(out + '.gen.gz', 'wt') as handle:
        handle.write('22 rs123 100 A C 1 0 0\n')
elif '--gen' in args or '--bfile' in args:
    out = value('--out')
    write(out + '.bed', 'bed\n')
    write(out + '.bim', '22 rs123 0 100 A C\n')
    write(out + '.fam', 'F I 0 0 1 -9\n')
else:
    raise SystemExit('unexpected mock arguments: ' + repr(args))
