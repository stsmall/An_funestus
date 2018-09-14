#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 12:43:47 2018

@author: scott
"""

from __future__ import print_function
from __future__ import division

import argparse
import allel
import numpy as np
import os
import h5py

parser = argparse.ArgumentParser()
parser.add_argument("--vcfFile", required=True,
                    help="vcfFile")
parser.add_argument("--cM", required=True,
                    help="physical position in cM")
parser.add_argument('-p', "--pops", nargs='+', action="append",
                    required=True, help="index of test pop A and B"
                    "-p 1 2 3 -p 6 7 8")
parser.add_argument('-o', "--outgroup", nargs='+', action="append",
                    help="index of outgroup")
args = parser.parse_args()


def loadvcf(vcfFile):
    """Reads VCF using scikit-allel, object is stored as pandas DF
    """
    print("loading vcf file...")
    print("using scitkit allele version:", allel.__version__)
    h5 = "{}h5".format(vcfFile.strip("vcf"))
    if os.path.isfile(h5):
        callset = h5py.File(h5, 'r')
    else:
        callset = allel.read_vcf(vcfFile)
        print("creating h5 for faster loading")
        allel.vcf_to_hdf5(vcfFile, h5)
    return(callset)


def filterGT(callset, outgroup):
    """Count patterns from VCF
    """
    gt = allel.GenotypeArray(callset['calldata/GT'])
    p = callset['variants/POS']
    pos = allel.SortedIndex(p)
    acs = gt[:, outgroup].count_alleles(max_allele=1)
    flt = acs.is_segregating()  # needs to be segregating in the outgroup
    gt = gt.compress(flt, axis=0)
    pos = pos[flt]
    return(gt, pos)


def countAlleles(gt, pops, outgrp):
    """
    """
    # make certain segregating in popC
    subpops = {"A": pops[0],
               "B": pops[1],
               "C": outgrp
               }
    acs = gt.count_alleles_subpops(subpops, max_allele=1)
    return(acs)


def make3PCLR(chrom, acs, cM, pos):
    """
    """
    f = open("{}.3pclrIn.txt".format(chrom), 'w')
    f.write("chr\tphypos\tgenpos\tmpopA\tnpopA\tmpopB\tnpopB\tmpopC\tnpopC\n")
    for p, c, a, b, c in zip(pos, cM, acs["A"], acs["B"], acs["C"]):
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, p, c, a[1],sum(a), b[1], sum(b), c[1], sum(c)))
    f.close()


if __name__ == "__main__":
    popset = args.pops
    pop_ix = [list(map(int, x)) for x in popset]
    outgroup_ix = [list(map(int, x)) for x in args.outgroup][0]
    callset = loadvcf(args.vcfFile)
    gt, pos = filterGT(callset, pop_ix, outgroup_ix)
    acs = countAlleles(gt, pop_ix, outgroup_ix)
