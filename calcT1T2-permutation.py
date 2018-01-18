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
parser.add_argument('-I', "--iterations", type=int, default=100,
                    help="number of iterations")
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
    callabledict = defaultdict(lambda: 0, d)
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
                        callabledict[chrom] += 1
                        for sample in range(9, len(x)):
                            ref = 0  # check for missing
                            alt = 0  # check for missing
                            gt = x[sample].split(":")[0]
                            ref += gt.count("0")
                            alt += gt.count("1")
                            if ref == 0 and alt == 0:
                                ref = -1
                                alt = -1
                            if "0/0" in polarize or "0|0" in polarize:
                                count_list.append([alt, ref])
                            elif "1/1" in polarize or "1|1" in polarize:
                                count_list.append([ref, alt])
                        if "0/1" not in polarize or "0|1" not in polarize or "1|0" not in polarize:
                            # if ancestral is not polymorphic
                            qdict[chrom][pos] = (count_list)
    return(qdict, q_ix, samplelist, callabledict)


def DfoilTble(t1t2dict, size, ntaxa):
    """
    """
    if ntaxa is 4:
        headers = ['AAAA', 'AABA', 'ABAA', 'ABBA',
                   'BAAA', 'BABA', 'BBAA', 'BBBA']
    elif ntaxa is 5:
        headers = ['AAAAA', 'AAABA', 'AABAA', 'AABBA',
                   'ABAAA', 'ABABA', 'ABBAA', 'ABBBA',
                   'BAAAA', 'BAABA', 'BABAA', 'BABBA',
                   'BBAAA', 'BBABA', 'BBBAA', 'BBBBA']
    d = open("dfoil.tbl", 'w')
    if size == 0:
        d.write("#chrom\tsites\t{}\n".format('\t'.join(headers)))
        for chrom in t1t2dict.keys():
            posdict = OrderedDict(sorted(t1t2dict[chrom].items()))
            divergence = []
            for pos in posdict.keys():
                posmean = np.mean(np.array(posdict[pos]), axis=0)
                divergence.append(posmean)
            div = np.array(divergence)
            sites = len(divergence)
            #sites = np.sum(div, axis=1)
            div_sum = np.sum(div, axis=0)
            divstr = map(str, div_sum)
            d.write("{}\t{}\t{}\n".format(chrom, sites, '\t'.join(divstr)))
    else:
        d.write("#chrom\tstart\tend\tsites\t{}\n".format('\t'.join(headers)))
        start = 1
        end = size
        for chrom in t1t2dict.keys():
            posdict = OrderedDict(sorted(t1t2dict[chrom].items()))
            divergence = []
            for pos in posdict.keys():
                if pos > end:
                    try:
                        # AAAA, AABA, ABAA, ABBA, BAAA, BABA, BBAA, BBBA
                        div = np.array(divergence)
                        sites = len(divergence)
                        #sites = np.sum(div, axis=1)
                        div_sum = np.sum(div, axis=0)
                        try:
                            divstr = map(str, div_sum)
                        except TypeError:
                            raise IndexError
                        d.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, sites, '\t'.join(divstr)))
                        divergence = []
                        start = end
                        end = end + size
                    except IndexError:
                        d.write("{}\t{}\t{}\t{}\t{}0\n".format(chrom, start, end, sites, '0\t'*7))
                        start = end
                        end = end + size
                else:
                    posmean = np.mean(np.array(posdict[pos]), axis=0)
                    divergence.append(posmean)
    d.close()
    return(None)


def foil4(vcfdict, quartet, q_ix, samplelist, iterations, callabledict):
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
        for its in range(iterations):
            i = np.random.choice(q_ix[0], 1)[0]
            j = np.random.choice(q_ix[1], 1)[0]
            k = np.random.choice(q_ix[2], 1)[0]
            n_AAAA = 0  #
            n_AABA = 0  #
            n_ABAA = 0
            n_ABBA = 0  #
            n_BAAA = 0
            n_BABA = 0  #
            n_BBAA = 0
            n_BBBA = 0  #
            countlist = []
            for pos in vcfdict[chrom].keys():
                marray = np.array(vcfdict[chrom][pos])
                m = np.array([marray[i-9], marray[j-9], marray[k-9], marray[-1]])
                if -1 not in m:
                    window = [0, 0, 0, 0, 0, 0, 0, 0]
                    header = ['AAAA', 'AABA', 'ABAA', 'ABBA', 'BAAA',
                              'BABA', 'BBAA', 'BBBA']
                    count = np.where(m == 0)
                    try:
                        count_sum = sum(count[1][0:3])  # only first 3
                    except IndexError:
                        import ipdb;ipdb.set_trace()
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
            if callabledict[chrom] > 0:
                # P1 P2 P3 O; BAAA, ABAA, BBAA
                t2_inner = (n_ABAA + n_BAAA) / 2
                t2_1 = t2_inner / callabledict[chrom]
                t1_1 = (t2_inner + n_BBAA) / callabledict[chrom]
                # t1se, t2se = blockSE(t1t2dict, 2, 4, 6)

                # P1 P3 P2 O; BAAA AABA BABA
                t2_inner = (n_BAAA + n_AABA) / 2
                t2_2 = t2_inner / callabledict[chrom]
                t1_2 = (t2_inner + n_BABA) / callabledict[chrom]
                # t1se, t2se = blockSE(t1t2dict, 4, 1, 5)

                # P2 P3 P1 O; ABAA AABA ABBA
                t2_inner = (n_ABAA + n_AABA) / 2
                t2_3 = t2_inner / callabledict[chrom]
                t1_3 = (t2_inner + n_ABBA) / callabledict[chrom]
                # t1se, t2se = blockSE(t1t2dict, 2, 1, 3)
            t1list.append([t1_1, t1_2, t1_3])
            t2list.append([t2_1, t2_2, t2_3])
            countlist.append([n_AAAA, n_AABA, n_ABAA, n_ABBA, n_BAAA,
                              n_BABA, n_BBAA, n_BBBA])
        # averages w/ SE
        reps = len(t1list)
        np.savetxt("t1array.out", t1list)
        np.savetxt("t2array.out", t2list)
        np.savetxt("counts.out", countlist)
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


def foil5(vcfdict, quartet, q_ix, samplelist, iterations):
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
        for its in range(iterations):
            i = np.random.choice(q_ix[0], 1)
            j = np.random.choice(q_ix[1], 1)
            k = np.random.choice(q_ix[2], 1)
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
    size = args.size
    qdict, q_ix, samplelist, calldict = loadvcf(vcfFile, quart, args.dlm)
    if len(quart) == 5:
        t1t2dict = foil5(qdict, quart, q_ix, samplelist, args.iterations, calldict)
    elif len(quart) == 4:
        t1t2dict = foil4(qdict, quart, q_ix, samplelist, args.iterations, calldict)
    else:
        raise ValueError("quartet must be 4 or 5 taxa")
    DfoilTble(t1t2dict, size, len(quart))
