#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 12:43:47 2018

maplength in cM
##contig=<ID=3L,length=46588628>
##contig=<ID=3R,length=43534138>
##contig=<ID=2L,length=44382598>
##contig=<ID=2R,length=54617495>
##contig=<ID=X,length=20138207>

@author: scott
"""

from __future__ import print_function
from __future__ import division

import argparse
import allel
import bisect
import os
import h5py
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--vcfFile", required=True,
                    help="vcfFile")
parser.add_argument("--cMMb", required=True,
                    help="rate between SNPs in cM/Mb")
parser.add_argument('-p', "--pops", nargs='+', action="append",
                    required=True, help="index of test pop A and B"
                    "-p 1 2 3 -p 6 7 8")
parser.add_argument('-o', "--outgroup", nargs='+', action="append",
                    help="index of outgroup")
parser.add_argument("--chrom", type=str, help="name of chrom for outfile")
parser.add_argument("--mapsize", type=int, help="mapsize", required=True)
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
    for p1, cM1, a1, b1, c1 in zip(pos, cM, acs["A"], acs["B"], acs["C"]):
        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, p1, cM1, a1[1],sum(a1), b1[1], sum(b1), c1[1], sum(c1)))
    f.close()


def makecMmap(cMFile, pos, size):
    """
    """
    snplist = []
    cMlist = []
    cMlist2 = []
    with open(cMFile, 'r') as cm:
        for line in cm:
            x = line.split()
            snplist.append(int(x[0]))
            cMlist.append(float(x[1]))
    for p1 in pos:
        cMlist2.append(np.interp(p1, snplist, cMlist))
#        if p1 <= snplist[0]:
#            cM = cMlist[0] * (p1 / snplist[0])
#        else:
#            ixr = bisect.bisect(snplist, p1)
#            ixl = ixr - 1
#            cM = cMlist[ixl] + cMlist[ixl] * (p1 / snplist[ixr])
##        if i == 0:
##            cM += cMMb * p1 / size
##        else:
##            cM += (cMMb * (p1 - pos[i-1])) / size
#        cMlist2.append(cM)
##            if (ixr - ixl) > 1:  # average between the SNPs
##                if i == 0:
##                    total_dist = p
##                else:
##                    total_dist =  p - pos[i-1]
##                rw = (snplist[ixr] - p) / total_dist
##                lw = (p - snplist[ixl]) / total_dist
##                cMMb = (cMMblist[ixl] * lw) + (cMMblist[ixr] * rw)
##            else:
#                cMMb = cMMblist[ixl]
    return(cMlist2)


if __name__ == "__main__":
    popset = args.pops
    pop_ix = [list(map(int, x)) for x in popset]
    outgroup_ix = [list(map(int, x)) for x in args.outgroup][0]
    callset = loadvcf(args.vcfFile)
    gt, pos = filterGT(callset, outgroup_ix)
    acs = countAlleles(gt, pop_ix, outgroup_ix)
    cM = makecMmap(args.cMMb, pos, args.mapsize)
    make3PCLR(args.chrom, acs, cM, pos)
