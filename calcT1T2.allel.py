#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 17:59:46 2018
calcT1T2.py -t trees -v vcf -g group --nodes
Calculates the T1 and T2 divergence times in a quartet
@author: stsmall
"""
from __future__ import print_function
from __future__ import division
import numpy as np
import allel
import h5py
from itertools import product
from collections import defaultdict
import os.path
import bisect
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, required=True,
                    help="vcf file of variants")
parser.add_argument('-g', "--groups", nargs='+', action="append",
                    required=True, help="index of group"
                    "-g 1 2 3 -g 6 7 8 -g 11 12 13")
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
    print("using scitkit allele version:", allel.__version__)
    h5 = "{}h5".format(vcfFile.strip("vcf"))
    if os.path.isfile(h5):
        callset = h5py.File(h5, 'r')
    else:
        callset = allel.read_vcf(vcfFile)
        print("creating h5 for faster loading")
        allel.vcf_to_hdf5(vcfFile, h5)
    # assert len(set(callset['variants/CHROM'])) == 1
    chrom = callset['variants/CHROM'][0]
    return(callset, chrom)


def calct1t2(AAAA, BAAA, ABAA, AABA, BBAA, ABBA, BABA, BBBA):
    """
    """
    t1list = []
    t2list = []
    sites = AAAA + BAAA + ABAA + AABA + BBAA + ABBA + BABA + BBBA
    if sites != 0:
        # BAAA,ABAA,BBAA
        t2_inner = (BAAA + ABAA) / 2
        t2 = t2_inner / sites
        t1 = (t2_inner + BBAA) / sites
        t1list.append(t1)
        t2list.append(t2)
        # ABAA, AABA, ABBA
        t2_inner = (ABAA + AABA) / 2
        t2 = t2_inner / sites
        t1 = (t2_inner + ABBA) / sites
        t1list.append(t1)
        t2list.append(t2)
        # BAAA,AABA, BABA
        t2_inner = (BAAA + AABA) / 2
        t2 = t2_inner / sites
        t1 = (t2_inner + BABA) / sites
        t1list.append(t1)
        t2list.append(t2)
    else:
        t1list = [0, 0, 0]
        t2list = [0, 0, 0]
    return(t1list, t2list)


def countPattern(callset, sample_ix, outgroup):
    """Count patterns for all samples
    """
    print("counting patterns in file...")
    gt = allel.GenotypeArray(callset['calldata/GT'])
    pos = allel.SortedIndex(callset['variants/POS'])
    # remove any sites where outgroup is ./. or 0/1
    keep = gt[:, outgroup].is_hom() & gt.count_alleles().is_biallelic()
    gt = gt.compress(keep, axis=0)
    pos = pos[keep]
    # permute among all sample indexes, list of lists
    # [[1,2,3,4,5],[6,7,8,9],[12,14,15,16]]
    t1t2dict = defaultdict(list)
    windict = {}
    permute = 1
    g1, g2, g3 = sample_ix
    quartet = list(product(g1, g2, g3))
    print("total number of combinations: {}".format(len(quartet)))
    for quart in quartet:
        print("permutation number {}".format(permute))
        i, j, k = quart
        gt_sub = gt.take([i, j, k, outgroup], axis=1)
        keep = gt_sub.is_hom().all(axis=1)
        gt_sub = gt_sub.compress(keep, axis=0)
        pos_sub = pos[keep]
        count_array = gt_sub.is_hom_alt()
        pattern_array = np.packbits(count_array, axis=1)
        calc_patterns = np.unique(pattern_array, return_counts=True)
        d = {n: calc_patterns[1][i] for i, n in enumerate(calc_patterns[0])}
        # total counts
        AAAA = d.get(0, 0) + d.get(240, 0)  # FFFF TTTT 240 and 0
        BAAA = d.get(112, 0) + d.get(128, 0)  # FTTT + TFFF 112 and 128
        ABAA = d.get(176, 0) + d.get(64, 0)  # TFTT + FTFF 176 and 64
        AABA = d.get(208, 0) + d.get(32, 0)  # TTFT + FFTF 208 and 32
        BBAA = d.get(48, 0) + d.get(192, 0)  # FFTT + TTFF 48 and 192
        ABBA = d.get(144, 0) + d.get(96, 0)  # TFFT + FTTF 144 and 96
        BABA = d.get(80, 0) + d.get(160, 0)  # FTFT + TFTF 80 and 160
        BBBA = d.get(224, 0) + d.get(16, 0)  # FFFT + TTTF 224 and 16
        # t1t2 calc
        t1, t2 = calct1t2(AAAA, BAAA, ABAA, AABA, BBAA, ABBA, BABA, BBBA)
        t1t2dict["t1"].append(t1)
        t1t2dict["t2"].append(t2)
        # windows
        windict[permute] = (pos_sub, pattern_array)
        permute += 1
    return(t1t2dict, windict)


def windowPattern(windict, size, chrom):
    """
    """
    print("printing patterns to file...")
    patterndict = defaultdict(list)
    # windowT1 = defaultdict(list)
    # windowT2 = defaultdict(list)
    for k in windict.keys():
        start = 0
        step = size
        end = start + size
        try:
            pos_array = windict[k][0]
            last = pos_array[-1]
        except IndexError:
            continue
        while end < last:
            start_ix = bisect.bisect_left(pos_array, start)
            end_ix = bisect.bisect_left(pos_array, end)
            try:
                pattern_array = windict[k][1][start_ix:end_ix]
                calc_patterns = np.unique(pattern_array, return_counts=True)
                d = {n: calc_patterns[1][i] for i, n in enumerate(calc_patterns[0])}
                # total counts
                AAAA = d.get(0, 0) + d.get(240, 0)  # FFFF TTTT 240 and 0
                BAAA = d.get(112, 0) + d.get(128, 0)  # FTTT + TFFF 112 and 128
                ABAA = d.get(176, 0) + d.get(64, 0)  # TFTT + FTFF 176 and 64
                AABA = d.get(208, 0) + d.get(32, 0)  # TTFT + FFTF 208 and 32
                BBAA = d.get(48, 0) + d.get(192, 0)  # FFTT + TTFF 48 and 192
                ABBA = d.get(144, 0) + d.get(96, 0)  # TFFT + FTTF 144 and 96
                BABA = d.get(80, 0) + d.get(160, 0)  # FTFT + TFTF 80 and 160
                BBBA = d.get(224, 0) + d.get(16, 0)  # FFFT + TTTF 224 and 16
                # t1, t2 = calct1t2(AAAA, BAAA, ABAA, AABA, BBAA, ABBA, BABA, BBBA)
            except IndexError:
                pattern_array = windict[k][1][start_ix:]
                calc_patterns = np.unique(pattern_array, return_counts=True)
                d = {n: calc_patterns[1][i] for i, n in enumerate(calc_patterns[0])}
                # total counts
                AAAA = d.get(0, 0) + d.get(240, 0)  # FFFF TTTT 240 and 0
                BAAA = d.get(112, 0) + d.get(128, 0)  # FTTT + TFFF 112 and 128
                ABAA = d.get(176, 0) + d.get(64, 0)  # TFTT + FTFF 176 and 64
                AABA = d.get(208, 0) + d.get(32, 0)  # TTFT + FFTF 208 and 32
                BBAA = d.get(48, 0) + d.get(192, 0)  # FFTT + TTFF 48 and 192
                ABBA = d.get(144, 0) + d.get(96, 0)  # TFFT + FTTF 144 and 96
                BABA = d.get(80, 0) + d.get(160, 0)  # FTFT + TFTF 80 and 160
                BBBA = d.get(224, 0) + d.get(16, 0)  # FFFT + TTTF 224 and 16
                # t1, t2 = calct1t2(AAAA, BAAA, ABAA, AABA, BBAA, ABBA, BABA, BBBA)
            # windowT1[end].append(t1)
            # windowT2[end].append(t2)
            patterndict[end].append((AAAA, BAAA, ABAA, AABA, BBAA, ABBA, BABA, BBBA))
            start += step
            end += step
    # write to file
    wfile = open("dfoil.tbl", 'w')
    headers = ['AAAA', 'AABA', 'ABAA', 'ABBA',
               'BAAA', 'BABA', 'BBAA', 'BBBA']
    # headers = ["ABC", "BCA", "ACB", "AB", "BC", "AC"]
    wfile.write("chrom\tpos\t{}\n".format('\t'.join(headers)))
    ordered_keys = sorted(list(patterndict.keys()))
    for pos in ordered_keys:
        count_mean = np.mean(list(zip(*patterndict[pos])), axis=1)
        wfile.write("{}\t{}\t{}\n".format(chrom, pos, "\t".join(map(str, count_mean))))
    wfile.close()
    return(None)


if __name__ == "__main__":
    groups = args.groups
    quart_ix = [list(map(int, x)) for x in groups]
    outgroup_ix = int(args.outgroup)
    vcfFile = args.vcfFile
    size = args.size
    callset, chrom = loadvcf(vcfFile)
    t1t2dict, windict = countPattern(callset, quart_ix, outgroup_ix)
    i, j, k = list(zip(*t1t2dict["t1"]))
    m, n, p = list(zip(*t1t2dict["t2"]))
    print("(AB)C: {}".format(np.mean(i)))
    print("(BC)A: {}".format(np.mean(j)))
    print("(AC)B: {}".format(np.mean(k)))
    print("AB: {}".format(np.mean(m)))
    print("BC: {}".format(np.mean(n)))
    print("AC: {}".format(np.mean(p)))
    windowPattern(windict, size, chrom)
