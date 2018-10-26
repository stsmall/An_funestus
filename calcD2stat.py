#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 16:48:06 2018

@author: scott
"""

from __future__ import print_function
from __future__ import division
import numpy as np
from ete3 import PhyloTree
from itertools import combinations
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', "--treefile", type=str, required=True,
                    help="treefile in newick, 1 per line")
parser.add_argument('-g', "--groups", nargs='+', action='append',
                    required=True, help="quartet of species to calculate,"
                    " assumes form: P1 P2 P3. can be given multiple times")
parser.add_argument('-w', "--windows", type=str,
                    help="coordinates for each tree")
parser.add_argument("--dlm", type=str, default="_",
                    help="delimeter denoting species")
parser.add_argument("-o", "--outgroup", type=str,
                    help="outgroup species for rooting")
args = parser.parse_args()


def LoadTrees(treefile, outgroup, dlm):
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
    with open(treefile, 'r') as newick:
        for line in newick:
            if not line.startswith("NA"):
                t = PhyloTree(line)
                t.set_species_naming_function(lambda node: node.name.split(dlm)[0])
                t.set_outgroup( t&outgroup )
                treelist.append(t)
    return(treelist)


def cMono(tree, taxon):
    """Checks if samples are monophyletic in tree
    """
    return(tree.check_monophyly(values=taxon, target_attr="species")[0])


def pairwiseDistance(tree, taxon):
    """Iterative distance for each individual regardless of monophyly
    """
    pwlist = []
    sp1 = tree.search_nodes(species=taxon[0])
    sp2 = tree.search_nodes(species=taxon[1])
    for s1 in sp1:
        for s2 in sp2:
            pwlist.append(tree.get_distance(s1, s2))
    return(np.nanmean(pwlist))


def AgeAndSupport(treelist, taxon):
    """Calculates the support and node age if groups in taxon are monophyletic
    """
    # only A/C when (A,B) and (B,C)
    A = [taxon[0]]
    B = [taxon[1]]
    C = [taxon[2]]
    AC_AB = []
    AC_BC = []
    AB_AB = []
    BC_BC = []
    for t in treelist:
        ass = t.search_nodes(species=A[0])
        bss = t.search_nodes(species=B[0])
        css = t.search_nodes(species=C[0])
        t.prune(ass+bss+css, preserve_branch_length=True)
        if cMono(t, A+B):
            AC_AB.append(pairwiseDistance(t, A+C))
            AB_AB.append(pairwiseDistance(t, A+B))
        else:
            AC_AB.append(0)
            AB_AB.append(0)
        if cMono(t, B+C):
            AC_BC.append(pairwiseDistance(t, A+C))
            BC_BC.append(pairwiseDistance(t, B+C))
        else:
            AC_BC.append(0)
            BC_BC.append(0)
    return(AC_AB, AC_BC, AB_AB, BC_BC)


if __name__ == "__main__":
    quarts = args.groups
    quart = quarts[0]
    taxon = [quart[0], quart[1], quart[2]]
    treelist = LoadTrees(args.treefile, args.outgroup, args.dlm)
    dac_ab, dac_bc, dab_ab, dbc_bc= AgeAndSupport(treelist, taxon)
    # calculate sliding window by 100 trees or such
    i = 0
    step = 100
    j = step
    f = open("D2_clust.txt", 'w')
    while j < len(dac_ab):
        d2 = np.mean(dac_ab[i:j]) - np.mean(dac_bc[i:j])
        f.write("{}\n".format(d2))
        i = j
        j += step
    f.close()
    print("D2: {}".format(np.mean(dac_ab) - np.mean(dac_bc)))
    i = 0
    step = 100
    j = step
    f = open("D1_clust.txt", 'w')
    while j < len(dab_ab):
        d1 = np.mean(dab_ab[i:j]) - np.mean(dbc_bc[i:j])
        f.write("{}\n".format(d1))
        i = j
        j += step
    f.close()
    print("D1: {}".format(np.mean(dab_ab) - np.mean(dbc_bc)))
