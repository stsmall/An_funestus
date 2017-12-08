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
import ipdb
import numpy as np
import argparse
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, required=True,
                    help="vcf file of variants")
parser.add_argument('-g', "--groups", nargs='+',
                    required=True, help="quartet of species to calculate,"
                    " assumes form: P1 P2 P3. can be given multiple times")
parser.add_argument("--dlm", type=str, default=".",
                    help="delimeter denoting species")
args = parser.parse_args()

# TODO: allele freq for D rather than count


def loadvcf(vcFile, quart, dlm):
    """Creates a dictionary object from a vcffile only including species in the
    given quartet.
    """
    print("loading vcf file...")
    qdict = defaultdict(dict)
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                samplelist = line.strip().split()
                q_ix = []
                for q in quart:
                    q_ix.append([i for i, x in enumerate(samplelist) if q == x.split(dlm)[0]])
            elif not line.startswith("##"):
                x = line.strip().split()
                if "," in x[4]:
                    # skip triallelic sites
                    continue
                else:
                    chrom = x[0]
                    pos = x[1]
                    count_list = []
                    polarize = x[q_ix[-1][0]].split(":")[0]
                    if "." not in polarize:
                        for sample in range(9, len(x)):
                            ref = 0  # check for missing
                            alt = 0  # check for missing
                            gt = x[sample].split(":")[0]
                            ref += gt.count("0")
                            alt += gt.count("1")
                            if ref == 0 and alt == 0:
                                ref = -1
                                alt = -1
                            if "0/0" in polarize:
                                count_list.append([alt, ref])
                            elif "1/1" in polarize:
                                count_list.append([ref, alt])
                        if "0/1" not in polarize:
                            # if ancestral is not polymorphic
                            qdict[chrom][pos] = (count_list)
    return(qdict, q_ix, samplelist)


def foil4(vcfdict, quartet, q_ix, samplelist):
    """Calculates the divergence between (1,2) as:
        T2 = (1/N) * ((n_ABAA + n_BAAA) / 2).
      Calculates the divergence between (1,2),3 as:
        T1 = (1/N) * (T2 + n_BBAA)

    Parameters
    ------
    vcfdict: dict, obj from loadvcf
    quartet: list, list of groups

    Returns
    ------
    t1t2dict: dict, chrom : pos : (counts)

    """
    print("calculating divergence times for quartet: {}...".format(quartet))
    p1, p2, p3, p4 = quartet
    t1t2dict = defaultdict(lambda: defaultdict(lambda: []))
    for chrom in vcfdict.keys():
        t1list = []
        t2list = []
        indslist = []
        for i in q_ix[0]:
            for j in q_ix[1]:
                for k in q_ix[2]:
                    n_AAAA = 0  #
                    n_AABA = 0  #
                    n_ABAA = 0
                    n_ABBA = 0  #
                    n_BAAA = 0
                    n_BABA = 0  #
                    n_BBAA = 0
                    n_BBBA = 0  #
                    callable_pos = 0
                    countlist = []
                    for pos in vcfdict[chrom].keys():
                        marray = np.array(vcfdict[chrom][pos])
                        m = np.array([marray[i-9], marray[j-9], marray[k-9], marray[-1]])
                        if -1 not in m:
                            window = [0, 0, 0, 0, 0, 0, 0, 0]
                            header = ['AAAA', 'AABA', 'ABAA', 'ABBA', 'BAAA',
                                      'BABA', 'BBAA', 'BBBA']
                            callable_pos += 1
                            count = np.where(m == 0)
                            count_sum = sum(count[1][0:3])  # only first 3
                            count_len = len(count[1])  # 4 zeros
                            if count_len == 4:
                                if m[3, 1] != 0:
                                    if count_sum == 0:
                                        n_AAAA += 1
                                        window[header.index('AAAA')] = 1
                                    elif count_sum == 1:
                                        iix = np.where(count[1] == 1)[0]
                                        if 2 in iix:
                                            n_AABA += 1
                                            window[header.index('AABA')] = 1
                                        elif 1 in iix:
                                            n_ABAA += 1
                                            window[header.index('ABAA')] = 1
                                        elif 0 in iix:
                                            n_BAAA += 1
                                            window[header.index('BAAA')] = 1
                                    elif count_sum == 2:
                                        # two zeros
                                        iix = np.where(count[1] == 1)[0]
                                        if 0 in iix and 2 in iix:
                                            n_BABA += 1
                                            window[header.index('BABA')] = 1
                                        elif 0 in iix and 1 in iix:
                                            n_BBAA += 1
                                            window[header.index('BBAA')] = 1
                                        elif 1 in iix and 2 in iix:
                                            n_ABBA += 1
                                            window[header.index('ABBA')] = 1
                                    elif count_sum == 3:
                                        n_BBBA += 1
                                        window[header.index('BBBA')] = 1
                                    else:
                                        raise ValueError("unknown pattern")
                                    t1t2dict[chrom][int(pos)].append(window)
                    # 'AAAA', 'AABA', 'ABAA', 'ABBA', 'BAAA', 'BABA', 'BBAA', 'BBBA'
                    # 0        1        2      3        4       5      6       7
                    if callable_pos > 0:
                        # P1 P2 P3 O; BAAA, ABAA, BBAA
                        t2_inner = (n_ABAA + n_BAAA) / 2
                        t2_1 = t2_inner / callable_pos
                        t1_1 = (t2_inner + n_BBAA) / callable_pos
                        # t1se, t2se = blockSE(t1t2dict, 2, 4, 6)

                        # P1 P3 P2 O; BAAA AABA BABA
                        t2_inner = (n_BAAA + n_AABA) / 2
                        t2_2 = t2_inner / callable_pos
                        t1_2 = (t2_inner + n_BABA) / callable_pos
                        # t1se, t2se = blockSE(t1t2dict, 4, 1, 5)

                        # P2 P3 P1 O; ABAA AABA ABBA
                        t2_inner = (n_ABAA + n_AABA) / 2
                        t2_3 = t2_inner / callable_pos
                        t1_3 = (t2_inner + n_ABBA) / callable_pos
                        # t1se, t2se = blockSE(t1t2dict, 2, 1, 3)
                    t1list.append([t1_1, t1_2, t1_3])
                    t2list.append([t2_1, t2_2, t2_3])
                    inds = (samplelist[marray[i-9]], samplelist[marray[j-9]],
                            samplelist[marray[k-9]], samplelist[marray[-1]])
                    indslist.append(inds)
                    countlist.append([n_AAAA, n_AABA, n_ABAA, n_ABBA, n_BAAA,
                                      n_BABA, n_BBAA, n_BBBA])
        # averages w/ SE
        reps = len(t1list)
        ipdb.set_trace()
        np.savetxt("t1array.out", t1list)
        np.savetxt("t2array.out", t2list)
        np.savetxt("indslist.out", indslist)
        np.savetxt("counts.out", countlist)
        ipdb.set_trace()
        t1_1, t1_2, t1_3 = zip(*t1list)
        t2_1, t2_2, t2_3 = zip(*t2list)
        t1se = (np.std(t1_1)) / np.sqrt(reps)
        t2se = (np.std(t2_1)) / np.sqrt(reps)
        t1ase = (np.std(t1_2)) / np.sqrt(reps)
        t2ase = (np.std(t2_2)) / np.sqrt(reps)
        t1bse = (np.std(t1_3)) / np.sqrt(reps)
        t2bse = (np.std(t2_3)) / np.sqrt(reps)
        countmean = np.mean(countlist, axis=0)
        print("{}\n{}\n".format(header, countmean))
        print("{}\t({},{}),{} : {} +-{}\t({},{}) : {} +-{}\n".format(
              chrom, p1, p2, p3, np.mean(t1_1), t1se, p1, p2, np.mean(t2_1), t2se))
        print("{}\t({},{}),{} : {}+-{}\t({},{}) : {}+-{}\n".format(
                chrom, p1, p3, p2, np.mean(t1_2), t1ase, p1, p3, np.mean(t2_2), t2ase))
        print("{}\t({},{}),{} : {}+-{}\t({},{}) : {}+-{}\n".format(
                chrom, p2, p3, p1, np.mean(t1_3), t1bse, p2, p3, np.mean(t2_3), t2bse))
    return(t1t2dict)


def foil5(vcfdict, quartet, q_ix, samplelist):
    """Count pattern in VCF for DFOIL

    Parameters
    ------
    vcfdict: dict, obj from loadvcf
    quartet: list, list of groups

    Returns
    ------
    t1t2dict: dict, chrom : pos : (counts)

    """
    print("calculating divergence times for quartet: {}...".format(quartet))
    t1t2dict = defaultdict(lambda: defaultdict(lambda: []))
    for chrom in vcfdict.keys():
        for i in q_ix[0]:
            for j in q_ix[1]:
                for k in q_ix[2]:
                    n_AAAAA = 0
                    n_AAABA = 0
                    n_AABAA = 0
                    n_AABBA = 0
                    n_ABAAA = 0
                    n_ABABA = 0
                    n_ABBAA = 0
                    n_ABBBA = 0
                    n_BAAAA = 0
                    n_BAABA = 0
                    n_BABAA = 0
                    n_BABBA = 0
                    n_BBAAA = 0
                    n_BBABA = 0
                    n_BBBAA = 0
                    n_BBBBA = 0
                    callable_pos = 0
                    for pos in vcfdict[chrom].keys():
                        marray = np.array(vcfdict[chrom][pos])
                        m = np.array([marray[i-9], marray[j-9], marray[k-9], marray[-1]])
                        if -1 not in m:
                            window = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                            header = ['AAAAA', 'AAABA', 'AABAA', 'AABBA',
                                      'ABAAA', 'ABABA', 'ABBAA', 'ABBBA',
                                      'BAAAA', 'BAABA', 'BABAA', 'BABBA',
                                      'BBAAA', 'BBABA', 'BBBAA', 'BBBBA']
                            callable_pos += 1
                            count = np.where(m == 0)
                            count_sum = sum(count[1][0:4])  # sum only first 3 entries
                            count_len = len(count[1])  # 4 zeros
                            if count_len == 5:
                                if m[4, 1] != 0:
                                    if count_sum == 0:
                                        # 'AAAAA'
                                        n_AAAAA += 1
                                        window[header.index('AAAAA')] = 1
                                    elif count_sum == 1:
                                        # AAABA, AABAA, ABAAA, BAAAA, BBBBA
                                        iix = np.where(count[1] == 1)[0]
                                        if 0 in iix:
                                            n_BAAAA += 1
                                            window[header.index('BAAAA')] = 1
                                        elif 1 in iix:
                                            n_ABAAA += 1
                                            window[header.index('ABAAA')] = 1
                                        elif 2 in iix:
                                            n_AABAA += 1
                                            window[header.index('AABAA')] = 1
                                        elif 3 in iix:
                                            n_AAABA += 1
                                            window[header.index('AAABA')] = 1
                                    elif count_sum == 2:
                                        # AABBA, ABABA, ABBAA, BAABA, BABAA, BBAAA
                                        iix = np.where(count[1] == 1)[0]
                                        if 2 in iix and 3 in iix:
                                            n_AABBA += 1
                                            window[header.index('AABBA')] = 1
                                        elif 1 in iix and 3 in iix:
                                            n_ABABA += 1
                                            window[header.index('ABABA')] = 1
                                        elif 1 in iix and 2 in iix:
                                            n_ABBAA += 1
                                            window[header.index('ABBAA')] = 1
                                        elif 0 in iix and 3 in iix:
                                            n_BAABA += 1
                                            window[header.index('BAABA')] = 1
                                        elif 0 in iix and 2 in iix:
                                            n_BABAA += 1
                                            window[header.index('BABAA')] = 1
                                        elif 0 in iix and 1 in iix:
                                            n_BBAAA += 1
                                            window[header.index('BBAAA')] = 1
                                    elif count_sum == 3:
                                        # ABBBA, BABBA, BBABA, BBBAA
                                        iix = np.where(count[1] == 1)[0]
                                        if 1 in iix and 2 in iix and 3 in iix:
                                            n_ABBBA += 1
                                            window[header.index('ABBBA')] = 1
                                        elif 0 in iix and 2 in iix and 3 in iix:
                                            n_BABBA += 1
                                            window[header.index('BABBA')] = 1
                                        elif 0 in iix and 1 in iix and 3 in iix:
                                            n_BBABA += 1
                                            window[header.index('BBABA')] = 1
                                        elif 0 in iix and 1 in iix and 2 in iix:
                                            n_BBBAA += 1
                                            window[header.index('BBBAA')] = 1
                                    elif count_sum == 4:
                                        n_BBBBA += 1
                                        window[header.index('BBBBA')] = 1
                                    else:
                                        raise ValueError("pattern not recognized")
                                    t1t2dict[chrom][int(pos)].append(window)
    return(t1t2dict)


if __name__ == "__main__":
    quart = args.groups
    vcfFile = args.vcfFile
    qdict, q_ix, samplelist = loadvcf(vcfFile, quart, args.dlm)
    if len(quart) == 5:
        t1t2dict = foil5(qdict, quart, q_ix, samplelist)
    elif len(quart) == 4:
        t1t2dict = foil4(qdict, quart, q_ix, samplelist)
    else:
        raise ValueError("quartet must be 4 or 5 taxa")