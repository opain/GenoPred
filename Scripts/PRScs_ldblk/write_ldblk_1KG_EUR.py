#!/usr/bin/env python

"""
Read LD blocks and write as hdf5

"""

import scipy as sp

import h5py
import argparse

p = argparse.ArgumentParser()
p.add_argument('-blkdir', required=True)
p.add_argument('-out', required=True)
args = p.parse_args()


BLK_DIR = args.blkdir
OUT_DIR = args.out


with open(BLK_DIR + '/LD_Blocks/blk_chr') as ff:
    blk_chr = [int(line.strip()) for line in ff]


with open(BLK_DIR + '/LD_Blocks/blk_size') as ff:
    blk_size = [int(line.strip()) for line in ff]


n_blk = len(blk_chr)

n_chr = max(blk_chr)


for chrom in range(1,n_chr+1):

    print('... parse chomosome %d ...' % chrom)
    chr_name = OUT_DIR + '/ldblk_1kg_chr' + str(chrom) + '.hdf5'
    hdf_chr = h5py.File(chr_name, 'w')
    blk_cnt = 0
    for blk in range(n_blk):
        if blk_chr[blk] == chrom:
            if blk_size[blk] > 0:
                blk_file = BLK_DIR + '/LD_Blocks/1KGPhase3.w_hm3.EUR.Block_' + str(blk+1) + '.ld'
                with open(blk_file) as ff:
                    ld = [[float(val) for val in (line.strip()).split()] for line in ff]
                print('blk %d size %s' % (blk+1, sp.shape(ld)))
                snp_file = BLK_DIR + '/LD_Blocks/Block_' + str(blk+1) + '.snplist'
                with open(snp_file) as ff:
                    snplist = [line.strip() for line in ff]
            else:
                ld = []; snplist = []
            blk_cnt += 1
            hdf_blk = hdf_chr.create_group('blk_%d' % blk_cnt)
            hdf_blk.create_dataset('ldblk', data=sp.array(ld), compression="gzip", compression_opts=9)
            hdf_blk.create_dataset('snplist', data=snplist, compression="gzip", compression_opts=9)

