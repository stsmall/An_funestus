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


def loadvcf(vcfile, quart):
    """Creates a dictionary object from a vcffile only including species in the
    given quartet.
    """
    return(vcfdict)


def calcT1(vcfdict, ABBB, BABB):
    """Calculates the divergence between (1,2),3.
    """
    t1 = ""
    return(t1)


def calcT2(vcfdict, quartet, size):
    """Calculates the divergence between (1,2).
    """
    t1list = []
    t2list = []

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
        raise ValueError("to calc node height need a tree file")
    vcfFile = args.vcffile
    quart = args.groups
    trees = args.treefile
    vcfdict = loadvcf(vcfFile)
    t1, t2 = calcT2(vcfdict, quart, args.size)
    treelist = loadtrees(trees)
    nh1, nh2 = nodeHeights(treelist, quart, args.windows)
