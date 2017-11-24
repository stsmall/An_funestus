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

import numpy as np
from ete3 import Tree
import argparse
from collections import defaultdict
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
parser.add_arguement("--nodes", action="store_true",
                     help="calculate node heights for a given quartet")
args = parser.parse_args()


def loadvcf(vcFile, quart):
    """Creates a dictionary object from a vcffile only including species in the
    given quartet.
    """
    qdict = defaultdict(dict)
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                sample = line.strip().split()
                q_ix = []
                for q in quart:
                    q_ix.append([i for i, x in enumerate(sample) if q in x])
            elif not line.startwith("##"):
                x = line.strip().split()
                chrom = x[0]
                pos = x[1]
                count_list = []
                for q in q_ix:
                    ref = 0
                    alt = 0
                    for s in q:
                        gt = x[s].split(":")[0]
                        ref += gt.count("0")
                        alt += gt.count("1")
                    count_list.append([ref, alt])
                qdict[chrom][pos] = (count_list)
    return(qdict)


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


    """
    t1list = []
    t2list = []
    for chrom in vcfdict.keys():
        n_ABAA = 0
        n_BAAA = 0
        n_BBAA = 0
        for pos in vcfdict[chrom].keys:


    return(t1list, t2list)


def loadtrees(treefile):
    """Loads trees from a treefile into an ete3 object
    """


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


if __name__ == "__main__":
    if args.nodes and not args.treefile:
        raise ValueError("to calc node heights need a tree file")
    vcfFile = args.vcffile
    quart = args.groups
    trees = args.treefile
    qdict = loadvcf(vcfFile, quart)
    t1, t2 = calcT2(qdict, quart, args.size)
    treelist = loadtrees(trees)
    nh1, nh2 = nodeHeights(treelist, quart, args.windows)
