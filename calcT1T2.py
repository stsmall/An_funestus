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
from ete3 import Tree
import argparse
from collections import defaultdict
from collections import OrderedDict
parser = argparse.ArgumentParser()
parser.add_argument('-t', "--treefile", type=str,
                    help="treefile in newick, 1 per line")
parser.add_argument('-v', "--vcffile", type=str, required=True,
                    help="vcf file of variants")
parser.add_argument('-g', "--groups", nargs='+',
                    help="quartet of species to calculate, assumes form: P1 P2"
                    "P3 O")
parser.add_argument('-s', "--size", type=int, default=0,
                    help="size of window for T1, T2 calculations")
parser.add_argument('-w', "--windows", type=str,
                    help="coordinates for each tree")
parser.add_argument("--nodes", action="store_true",
                    help="calculate node heights for a given quartet")
args = parser.parse_args()


def loadvcf(vcFile, quart):
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
                    q_ix.append([i for i, x in enumerate(sample) if q in x])
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


def t1t2slidingwindow(t1t2dict, size):
    """
    """
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
                    div = np.array(divergence)
                    sites = len(divergence)
                    # calc t2
                    t2_inner = sum(np.sum(div, axis=0)[0:2]) / 2
                    t2 = t2_inner / sites
                    # calc t1
                    t1 = (t2_inner + np.sum(div, axis=0)[2]) / sites
                    mid = (end - start) / 2
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start,
                                                              end, mid, t1,
                                                              t2))
                    divergence = []
                    start = end
                    end = end + size
                except IndexError:
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start,
                                                              end, mid, 0, 0))
                    start = end
                    end = end + size
            else:
                divergence.append(posdict[pos])
    f.close()
    return(None)


def calcT2(vcfdict, quartet, size):
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
        n_ABAA = 0
        n_BAAA = 0
        n_BBAA = 0
        callable_pos = 0
        for pos in vcfdict[chrom].keys():
            m = np.array(vcfdict[chrom][pos])
            if -1 not in m:
                window = [0, 0, 0]
                callable_pos += 1
                count_anc = np.sum(m, axis=0)[0]
                count_der = np.sum(m, axis=0)[1]
                if ((m[0, 0] == count_anc) and (m[0, 1] == 0)) or ((m[0, 1] == count_der) and (m[0, 0] == 0)):
                    n_BAAA += 1
                    window[0] = 1
                elif ((m[1, 0] == count_anc) and (m[1, 1] == 0)) or ((m[1, 1] == count_der) and (m[1, 0] == 0)):
                    n_ABAA += 1
                    window[1] = 1
                elif ((sum(m[0:2, 0]) == count_anc) and (sum(m[0:2, 1]) == 0)) or ((sum(m[0:2, 1]) == count_der) and (sum(m[0:2, 0]) == 0)):
                    n_BBAA += 1
                    window[2] = 1
                else:
                    pass
                t1t2dict[chrom][int(pos)] = tuple(window)
        if callable_pos > 0:
            t2_inner = (n_ABAA + n_BAAA) / 2
            t2 = t2_inner / callable_pos
            t1 = (t2_inner + n_BBAA) / callable_pos
            print("{}\t({},{}),{} : {}\t({},{}) : {}".format(chrom, p1, p2, p3,
                                                             t1, p1, p2, t2))
            t1dict[chrom].append(t1)
            t2dict[chrom].append(t2)
    if size != 0:
        t1t2slidingwindow(t1t2dict, size)
    return(t1dict, t2dict)


def loadtrees(treefile, outgroup):
    """Reads and stores phylogenetic trees from a file

    Parameters
    ------
    treefile: file, file of newick trees, 1 per line
    outgroup: str, last entry from quartet

    Returns
    ------
    treelist: obj, ete3 object of trees

    """
    print("loading trees...")
    treelist = []
    with open(treefile, 'r') as t:
        for line in t:
            if not line.startswith("NA"):
                t1 = Tree(line)
                t1.set_outgroup(outgroup)
                treelist.append(t1)
    return(treelist)


def nodeHeights(treedict, quart, windows):
    """Calculates the node heights in a set of trees.

    Parameters
    ------
    trees: ete3 object, returned from function loadtrees
    quart: list, list of topologies to calculate node heights
    windows: list, genome coordinates for which the trees were made

    Returns
    ------
    T1_T2: file, writes values of T1 and T2 for the topologies to a file for
           each tree.
    t1_avg: float, average value of T1 for topology
    t2_avg: float, average value of T2 for topology
    """
    print("calculating node heights for quartets: {}".format(quart))


if __name__ == "__main__":
    if args.nodes and not args.treefile:
        raise ValueError("to calc node heights need a tree file")
    vcfFile = args.vcffile
    quart = args.groups
    qdict = loadvcf(vcfFile, quart)
    t1, t2 = calcT2(qdict, quart, args.size)
    if args.treefile:
        treelist = loadtrees(args.treefile, quart[-1])
        if args.nodes:
            nh1, nh2 = nodeHeights(treelist, quart, args.windows)
