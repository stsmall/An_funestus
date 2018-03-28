#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 17:59:46 2018
calcT1T2.py -t trees -v vcf -g group --nodes
Calculates the T1 and T2 divergence times in a quartet
@author: stsmall
"""
import numpy as np
import allel
from itertools import combinations
from collections import defaultdict
from collections import OrderedDict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, required=True,
                    help="vcf file of variants")
parser.add_argument('-g', "--groups", nargs='+',
                    required=True, help="index of group"
                    "-g 1,2,3,4,5 -g 6,7,8,9,10 -g 11,12,13,14,15")
parser.add_argument('-o', "--outgroup", required=True,
                    help="index of outgroup")
parser.add_argument('-s', "--size", type=int, default=10000,
                    help="size of window for calculations")
parser.add_argument('-P', "--permutations", type=int, default=0,
                    help="number of permutations, default is all")
args = parser.parse_args()


def loadvcf(vcfFile):
    """Reads VCF using scikit-allel, object is stored as pandas DF
    """
    print("loading vcf file...")
    print(allel.__version__)
    callset = allel.read_vcf(vcfFile)
    # assert len(set(callset['variants/CHROM'])) == 1
    return(callset)


def windowPattern(pos, count_array, size):
    """
    """
    while ...
    count_array[pos.index(1), pos.index(10000)]
    AAAA = TTTT + FFFF
    BAAA = FTTT + TFFF
    ABAA = TFTT + FTFF
    AABA = TTFT + FFTF
    BBAA = FFTT + TTFF
    ABBA = TFFT + FTTF
    BABA = FTFT + TFTF
    BBBA = FFFT + TTTF


    return(None)


def countPattern(callset, sample_ix, outgroup, permutations, size):
    """Count patterns for all samples
    """
    gt = allel.GenotypeArray(callset['calldata/GT'])
    pos = gt.allel.sortedindex()
    # remove any sites where outgroup is ./. or 0/1
    outgroup_het = ~gt.take(outgroup).is_het()
    outgroup_miss = ~gt.is_missing()
    outgroup_mask = outgroup_het + outgroup_miss
    gt = gt.subset(outgroup_mask)
    pos = pos[outgroup_mask]
    # permute among all sample indexes, list of lists
    # [[1,2,3,4,5],[6,7,8,9],[12,14,15,16]]
    for i, j, k in sample_ix:
        gt2 = gt.take([i, j, k, outgroup])
        het_mask = ~gt2.is_het()  # remove het sites in advance
        miss_mask = ~gt2.is_missing()  # remove sites with missing
        count_mask = het_mask + miss_mask
        gt2 = gt2.subset(count_mask)
        pos2 = pos[count_mask]
        count_array = gt2.is_hom_alt()
        # total counts
        AAAA = TTTT + FFFF
        BAAA = FTTT + TFFF
        ABAA = TFTT + FTFF
        AABA = TTFT + FFTF
        BBAA = FFTT + TTFF
        ABBA = TFFT + FTTF
        BABA = FTFT + TFTF
        BBBA = FFFT + TTTF
        # windowed counts
        winarray = windowPattern(pos2, count_array, size)


def DfoilTble(t1t2dict, size, ntaxa):
    """
    """
   headers = ['AAAA', 'AABA', 'ABAA', 'ABBA',
              'BAAA', 'BABA', 'BBAA', 'BBBA']
   d.write("chrom\tstart\tend\tsites\t{}\n".format('\t'.join(headers)))

    ix = [[2, 4, 6], [1, 2, 3], [1, 4, 5]]
    for pat in ix:
        i, j, k = pat
        t2_inner = (divsum[i] + divsum[j]) / 2
        t2 = t2_inner / sites
        t1 = (t2_inner + divsum[k]) / sites
        print("{} {} {} : {}, {}".format(header[i], header[j], header[k], t1, t2))
    return(t1t2dict)




if __name__ == "__main__":
    quart = args.groups
    vcfFile = args.vcfFile
    size = args.size
    qdict, q_ix, samplelist, calldict = loadvcf(vcfFile, quart, args.dlm)
    if len(quart) == 5:
        t1t2dict = foil5(qdict, quart, q_ix, samplelist, args.iterations, calldict)
    elif len(quart) == 4:
        t1t2dict = foil4(qdict, quart, q_ix, samplelist, args.iterations, calldict)
    else:
        raise ValueError("quartet must be 4 or 5 taxa")
    DfoilTble(t1t2dict, size, len(quart))
