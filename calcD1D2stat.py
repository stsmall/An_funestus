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
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', "--treefile", type=str, required=True,
                    help="treefile in newick, 1 per line")
parser.add_argument('-g', "--groups", nargs='+', action='append',
                    required=True, help="quartet of species to calculate,"
                    " assumes form: P1 P2 P3. can be given multiple times")
parser.add_argument("--dlm", type=str, default="_",
                    help="delimeter denoting species")
parser.add_argument("-o", "--outgroup", type=str,
                    help="outgroup species for rooting if tree are unrooted")
parser.add_argument("--windows", type=int, default=1000, help="sliding windows")
parser.add_argument("--mono", action="store_true", help="enforce monophyly for" 
                    " trees with >1 individual per species")
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
                if outgroup:
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


def DistABC(treelist, taxon, mono):
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
        distBC = (pairwiseDistance(t, B+C))
        distAC = (pairwiseDistance(t, A+C))
        distAB = (pairwiseDistance(t, A+B))
        distAC = (pairwiseDistance(t, A+C))       
        if mono:
            if cMono(t, B+C):              
                if distBC < distAB and distBC < distAC:
                    AB_AB.append(distBC)
                    AC_AB.append(distAC)
            if cMono(t, A+B):
                if distAB < distBC and distAB < distAC:
                    AC_BC.append(distAC)
                    BC_BC.append(distAB)            
        else:
            if distBC < distAB and distBC < distAC:
                AB_AB.append(distBC)
                AC_AB.append(distAC)
            if distAB < distBC and distAB < distAC:
                AC_BC.append(distAC)
                BC_BC.append(distAB)            
    return(AC_AB, AC_BC, AB_AB, BC_BC)


if __name__ == "__main__":
    quarts = args.groups
    quart = quarts[0]
    taxon = [quart[0], quart[1], quart[2]]
    treelist = LoadTrees(args.treefile, args.outgroup, args.dlm)
    dac_ab, dac_bc, dab_ab, dbc_bc= DistABC(treelist, taxon, args.mono)
    A = quart[0][0]
    B = quart[1][0]
    C = quart[2][0]
    step = args.windows 
    
    # D1
    d1_m = np.mean(dab_ab) - np.mean(dbc_bc)
    # windowed D1
    d1SE = []
    # D1 sliding windows
    i = 0
    j = step
    f = open("{}{}{}.D1.{}.txt".format(A, B, C, step), 'w')
    while j < len(dab_ab):
        d1_win = np.mean(dab_ab[i:j]) - np.mean(dbc_bc[i:j])
        d1SE.append(d1_win)
        f.write("{}\n".format(d1_win))
        i = j
        j += step
    f.close()
    # D1 SE
    n = len(d1SE)
    try:
        sv = ((n - 1) / n) * np.sum((np.array(d1SE) - d1_m) ** 2)
    except ZeroDivisionError:
        se = 0.0000000001
    se = np.sqrt(sv)
    # print D1
    print("D1 ns from 0: speciation + introgression\nD1 sig +pos speciation followed by introgression\nincreasing D1 is more recent introgression")
    print("D1: {}, {}".format(d1_m, se))

    # D2
    d2_m = np.mean(dac_ab) - np.mean(dac_bc)
    # windowed D2
    d2SE = []
    # D2 calculate sliding window by 100 trees or such 
    i = 0
    j = step
    f = open("{}{}{}.D2.{}.txt".format(A, B, C, step), 'w')
    while j < len(dac_ab):
        d2_win = np.mean(dac_ab[i:j]) - np.mean(dac_bc[i:j])
        d2SE.append(d2_win)
        f.write("{}\n".format(d2_win))
        i = j
        j += step
    f.close()
    # D2 SE
    n = len(d2SE)
    try:
        sv = ((n - 1) / n) * np.sum((d2SE - d2_m) ** 2)
    except ZeroDivisionError:
        se = 0.0000000001
    se = np.sqrt(sv)
    # print D2
    print("D2 ns from 0: C->B\nD2 sig +pos B->C\nincreasing +pos w/ dist from speciation")
    print("D2: {}, {}".format(d2_m, se))
