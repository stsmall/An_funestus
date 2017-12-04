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
                sample = line.strip().split()
                q_ix = []
                for q in quart:
                    q_ix.append([i for i, x in enumerate(sample) if q == x.split(dlm)[0]])
                # randomly subsample q_ix to use only 1 individual
                q_ix_ind = [np.random.choice(i, 1) for i in q_ix]
                samplelist = [sample[i[0]] for i in q_ix_ind]
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
                        for q in q_ix_ind:
                            ref = 0  # check for missing
                            alt = 0  # check for missing
                            for s in q:
                                gt = x[s].split(":")[0]
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
                            qdict[chrom][pos] = (count_list)
    print("{}:{}\n{}:{}\n{}:{}\n{}:{}\n".format(quart[0], samplelist[0],
                                                quart[1], samplelist[1],
                                                quart[3], samplelist[3],
                                                quart[4], samplelist[4]))
    return(qdict)


def blockSE(t1t2dict, iix1, iix2, iix3, size=0, reps=10):
    """
    """
    t1sedict = {}
    t2sedict = {}
    if size == 0:
        for chrom in t1t2dict.keys():
            t1list = []
            t2list = []
            posdict = OrderedDict(sorted(t1t2dict[chrom].items()))
            sites = len(posdict.keys())
            for i in range(reps):
                pos = np.random.choice(posdict.keys(), sites, replace=True)
                divergence = []
                for p in pos:
                    divergence.append(posdict[p])
                # calc t1, t2
                div = np.array(divergence)
                div_sum = np.sum(div, axis=0)
                t2_inner = (div_sum[iix1] + div_sum[iix2]) / 2
                t2 = t2_inner / sites
                t1 = (t2_inner + div_sum[iix3]) / sites
                t1list.append(t1)
                t2list.append(t2)
            t1sedict[chrom] = (np.std(t1list)) / np.sqrt(reps)
            t2sedict[chrom] = (np.std(t2list)) / np.sqrt(reps)
    else:
        pass
        # TODO: block resampling
        # bin t1t2dict.keys() into size for block resampling accounting for ld
    return(t1sedict[chrom], t2sedict[chrom])


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
                divergence.append(posdict[pos])
            div = np.array(divergence)
            sites = len(divergence)
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
                        div_sum = np.sum(div, axis=0)
                        divstr = map(str, div_sum)
                        d.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, sites, '\t'.join(divstr)))
                        divergence = []
                        start = end
                        end = end + size
                    except IndexError:
                        d.write("{}\t{}\t{}\t{}\t{}0\n".format(chrom, start, end, sites, '0\t'*15))
                        start = end
                        end = end + size
                else:
                    divergence.append(posdict[pos])
    d.close
    return(None)


def foil4(vcfdict, quartet):
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
    t1t2dict = defaultdict(dict)
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
                header = ['AAAA', 'AABA', 'ABAA', 'ABBA', 'BAAA', 'BABA',
                          'BBAA', 'BBBA']
                callable_pos += 1
                count = np.where(m == 0)
                try:
                    count_sum = sum(count[1][0:3])  # sum only first 3 entries
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
                            raise ValueError("pattern not recognized")
                        t1t2dict[chrom][int(pos)] = tuple(window)
        #'AAAA', 'AABA', 'ABAA', 'ABBA', 'BAAA', 'BABA', 'BBAA', 'BBBA'
        # 0        1        2      3        4       5      6       7
        if callable_pos > 0:
            # P1 P2 P3 O; BAAA, ABAA, BBAA
            t2_inner = (n_ABAA + n_BAAA) / 2
            t2 = t2_inner / callable_pos
            t1 = (t2_inner + n_BBAA) / callable_pos
            t1se, t2se = blockSE(t1t2dict, 2, 4, 6)
            print("BAAA:{}\tABAA:{}\tBBAA:{}\tN:{}".format(n_BAAA, n_ABAA,
                                                           n_BBAA,
                                                           callable_pos))
            print("{}\t({},{}),{} : {}+-{}\t({},{}) : {}+-{}\n".format(chrom, p1, p2,
                                                               p3, t1, t1se,
                                                               p1, p2, t2,
                                                               t2se))
            # P1 P3 P2 O; BAAA AABA BABA
            t2_inner = (n_BAAA + n_AABA) / 2
            t2a = t2_inner / callable_pos
            t1a = (t2_inner + n_BABA) / callable_pos
            t1se, t2se = blockSE(t1t2dict, 4, 1, 5)
            print("BAAA:{}\tABAA:{}\tBBAA:{}\tN:{}".format(n_BAAA, n_AABA,
                                                           n_BABA,
                                                           callable_pos))
            print("{}\t({},{}),{} : {}+-{}\t({},{}) : {}+-{}\n".format(chrom, p1, p3,
                                                               p2, t1a, t1se, p1, p3,
                                                               t2a, t2se))
            # P2 P3 P1 O; ABAA AABA ABBA
            t2_inner = (n_ABAA + n_AABA) / 2
            t2b = t2_inner / callable_pos
            t1b = (t2_inner + n_ABBA) / callable_pos
            t1se, t2se = blockSE(t1t2dict, 2, 1, 3)
            print("BAAA:{}\tABAA:{}\tBBAA:{}\tN:{}".format(n_ABAA, n_AABA,
                                                           n_ABBA,
                                                           callable_pos))
            print("{}\t({},{}),{} : {}+-{}\t({},{}) : {}+-{}\n".format(chrom,
                                                                       p2, p3,
                                                                       p1, t1b,
                                                                       t1se, p2, p3,
                                                                       t2b, t2se))
    return(t1t2dict)


def foil5(vcfdict, quartet):
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
    t1t2dict = defaultdict(dict)
    for chrom in vcfdict.keys():
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
            m = np.array(vcfdict[chrom][pos])
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
                        t1t2dict[chrom][int(pos)] = tuple(window)
    return(t1t2dict)


if __name__ == "__main__":
    quart = args.groups
    vcfFile = args.vcfFile
    qdict = loadvcf(vcfFile, quart, args.dlm)
    if len(quart) == 5:
        t1t2dict = foil5(qdict, quart)
    elif len(quart) == 4:
        t1t2dict = foil4(qdict, quart)
    else:
        raise ValueError("quartet must be 4 or 5 taxa")
    DfoilTble(t1t2dict, args.size, len(quart))
