#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 15:24:22 2017
calcT1T2.py -t trees -v vcf -g group --nodes
Calculates the T1 and T2 divergence times in a quartet
@author: stsmall
"""
from __future__ import print_function
from __future__ import division
# from IPython.display import HTML
import numpy as np
import argparse
from collections import defaultdict
from collections import OrderedDict
parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, required=True,
                    help="vcf file of variants")
parser.add_argument('-g', "--groups", nargs='+',
                    required=True, help="quartet of species to calculate,"
                    " assumes form: P1 P2 P3. can be given multiple times")
parser.add_argument('-s', "--size", type=int, default=0,
                    help="size of window for T1, T2 calculations")
parser.add_argument("--dlm", type=str, default=".",
                    help="delimeter denoting species")
parser.add_argument("--dfoil", action="store_true",
                    help="returns table for dfoil")
args = parser.parse_args()


def loadvcf(vcFile, quart, dlm):
    """Creates a dictionary object from a vcffile only including species in the
    given quartet.
    """
    print("loading vcf file...")
    qdict = defaultdict(dict)
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                sample = line.strip().split()
                q_ix = []
                for q in quart:
                    q_ix.append([i for i, x in enumerate(sample) if q == x.split(dlm)[0]])
            elif not line.startswith("##"):
                x = line.strip().split()
                chrom = x[0]
                pos = x[1]
                count_list = []
                for q in q_ix:
                    ref = 0  # check for missing
                    alt = 0  # check for missing
                    for s in q:
                        gt = x[s].split(":")[0]
                        ref += gt.count("0")
                        alt += gt.count("1")
                    if ref == 0 and alt == 0:
                        ref = -1
                        alt = -1
                    count_list.append([ref, alt])
                qdict[chrom][pos] = (count_list)
    return(qdict)


def t1t2slidingwindow(t1t2dict, size, dfoil, ntaxa):
    """
    """
    if dfoil:
        ntaxa = 4
        if ntaxa is 4:
            headers = ['AAAA', 'AABA', 'ABAA', 'ABBA',
                       'BAAA', 'BABA', 'BBAA', 'BBBA']
        elif ntaxa is 5:
            headers = ['AAAAA', 'AAABA', 'AABAA', 'AABBA',
                       'ABAAA', 'ABABA', 'ABBAA', 'ABBBA',
                       'BAAAA', 'BAABA', 'BABAA', 'BABBA',
                       'BBAAA', 'BBABA', 'BBBAA', 'BBBBA']
        d = open("dfoil.tbl", 'w')
        d.write("#chrom\tposition\t{}\n".format('\t'.join(headers)))
    f = open("t1t2windowed.out", 'w')
    start = 1
    end = size
    f.write("chrom\tstart\tend\tmid\tt1\tt2\n")
    for chrom in t1t2dict.keys():
        posdict = OrderedDict(sorted(t1t2dict[chrom].items()))
        divergence = []
        for pos in posdict.keys():
            if pos > end:
                try:
                    # AAAA, AABA, ABAA, ABBA, BAAA, BABA, BBAA, BBBA
                    div = np.array(divergence)
                    sites = len(divergence)
                    div_sum = np.sum(div)
                    # calc t2
                    t2_inner = (div_sum[2] + div_sum[4]) / 2
                    t2 = t2_inner / sites
                    # calc t1
                    t1 = (t2_inner + div_sum[6]) / sites
                    mid = (end - start) / 2
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start,
                                                              end, mid, t1,
                                                              t2))
                    if dfoil:
                        d.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, mid, '\t'.join(div_sum)))
                    divergence = []
                    start = end
                    end = end + size
                except IndexError:
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start,
                                                              end, mid, 0, 0))
                    if dfoil:
                        d.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, mid, '\t'.join(div_sum)))
                    start = end
                    end = end + size
            else:
                divergence.append(posdict[pos])
    f.close()
    return(None)


def calcT1T2(vcfdict, quartet, size, dfoil):
    """Calculates the divergence between (1,2) as:
        T2 = (1/N) * ((n_ABAA + n_BAAA) / 2).
      Calculates the divergence between (1,2),3 as:
        T1 = (1/N) * (T2 + n_BBAA)

    Parameters
    ------
    vcfdict: dict, obj from loadvcf
    quartet: list, list of groups
    size: int, sliding window size

    Returns
    ------
    t1dict: dict, chrom : pos : t1
    t2dict: dict, chrom : pos : t2

    """
    print("calculating divergence times for quartet: {}...".format(quartet))
    p1, p2, p3, p4 = quartet
    t1t2dict = defaultdict(dict)
    t1dict = defaultdict(list)
    t2dict = defaultdict(list)
    for chrom in vcfdict.keys():
        n_AAAA = 0  #
        n_AABA = 0  #
        n_ABAA = 0
        n_ABBA = 0  #
        n_BAAA = 0
        n_BABA = 0  #
        n_BBAA = 0
        n_BBBA = 0  #
        callable_pos = 0
        for pos in vcfdict[chrom].keys():
            m = np.array(vcfdict[chrom][pos])
            if -1 not in m:
                window = [0, 0, 0, 0, 0, 0, 0, 0]
                # AAAA, AABA, ABAA, ABBA, BAAA, BABA, BBAA, BBBA
                callable_pos += 1
                count = np.where(m == 0)
                count_sum = sum(count[1])
                count_len = len(count[1])
                if count_len == 4:
                    if (count_sum == 0 or count_sum == 4):
                        n_AAAA += 1
                        window[0] = 1
                    elif (count_sum == 1 or count_sum == 3):
                        if count_sum == 1:
                            # find the 1
                            iix = np.where(count[1] == 1)[0]
                            if 2 in iix:
                                n_AABA += 1
                                window[1] = 1
                            elif 1 in iix:
                                n_ABAA += 1
                                window[2] = 1
                            elif 0 in iix:
                                n_BAAA += 1
                                window[4] = 1
                            else:
                                n_BBBA += 1
                                window[7] = 1
                        elif count_sum == 3:
                            # find the 0
                            iix = np.where(count[1] == 0)[0]
                            if 2 in iix:
                                n_AABA += 1
                                window[1] = 1
                            elif 1 in iix:
                                n_ABAA += 1
                                window[2] = 1
                            elif 0 in iix:
                                n_BAAA += 1
                                window[4] = 1
                            elif 3 in iix:
                                n_BBBA += 1
                                window[7] = 1
                    elif count_sum == 2:
                        # two zeros
                        iix = np.where(count[1] == 0)[0]
                        if 0 in iix and 2 in iix:
                            n_BABA += 1
                            window[5] = 1
                        elif 0 in iix and 1 in iix:
                            n_BBAA += 1
                            window[6] = 1
                        elif 1 in iix and 2 in iix:
                            n_ABBA += 1
                            window[3] = 1
                    else:
                        raise ValueError("pattern not recognized")
                        # import ipdb;ipdb.set_trace()
                t1t2dict[chrom][int(pos)] = tuple(window)

#        print("AAAA:{}".format(n_AAAA))
#        print("AABA:{}".format(n_AABA))
#        print("ABAA:{}".format(n_ABAA))
#        print("ABBA:{}".format(n_ABBA))
#        print("BAAA:{}".format(n_BAAA))
#        print("BABA:{}".format(n_BABA))
#        print("BBAA:{}".format(n_BBAA))
#        print("BBBA:{}".format(n_BBBA))

        if callable_pos > 0:
            # P1 P2 P3 O; BAAA, ABAA, BBAA
            t2_inner = (n_ABAA + n_BAAA) / 2
            t2 = t2_inner / callable_pos
            t1 = (t2_inner + n_BBAA) / callable_pos
            print("BAAA:{}\tABAA:{}\tBBAA:{}\tN:{}".format(n_BAAA, n_ABAA,
                                                           n_BBAA,
                                                           callable_pos))
            print("{}\t({},{}),{} : {}\t({},{}) : {}\n".format(chrom, p1, p2,
                                                               p3, t1, p1, p2,
                                                               t2))
            # P1 P3 P2 O; BAAA AABA BABA
            t2_inner = (n_BAAA + n_AABA) / 2
            t2a = t2_inner / callable_pos
            t1a = (t2_inner + n_BABA) / callable_pos
            print("BAAA:{}\tABAA:{}\tBBAA:{}\tN:{}".format(n_BAAA, n_AABA,
                                                           n_BABA,
                                                           callable_pos))
            print("{}\t({},{}),{} : {}\t({},{}) : {}\n".format(chrom, p1, p3,
                                                               p2, t1a, p1, p3,
                                                               t2a))
            # P2 P3 P1 O; ABAA AABA ABBA
            t2_inner = (n_ABAA + n_AABA) / 2
            t2b = t2_inner / callable_pos
            t1b = (t2_inner + n_ABBA) / callable_pos
            print("BAAA:{}\tABAA:{}\tBBAA:{}\tN:{}".format(n_ABAA, n_AABA,
                                                           n_ABBA,
                                                           callable_pos))
            print("{}\t({},{}),{} : {}\t({},{}) : {}\n".format(chrom, p2, p3,
                                                               p1, t1b, p2, p3,
                                                               t2b))
            t1dict[chrom].append(t1)
            t2dict[chrom].append(t2)
    if size != 0:
        t1t2slidingwindow(t1t2dict, size, dfoil, len(quartet))
    return(t1dict, t2dict)


if __name__ == "__main__":
    quart = args.groups
    vcfFile = args.vcfFile
    qdict = loadvcf(vcfFile, quart, args.dlm)
    t1, t2 = calcT1T2(qdict, quart, args.size, args.dfoil)
