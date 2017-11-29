#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 13:57:52 2017

@author: stsmall
"""
from __future__ import print_function
from __future__ import division
# from IPython.display import HTML
import numpy as np
from ete3 import PhyloTree
import argparse
from itertools import combinations
parser = argparse.ArgumentParser()
parser.add_argument('-t', "--treefile", type=str, required=True,
                    help="treefile in newick, 1 per line")
parser.add_argument('-g', "--groups", nargs='+',
                    required=True, help="quartet of species to calculate,"
                    " assumes form: P1 P2 P3. can be given multiple times")
parser.add_argument('-w', "--windows", type=str,
                    help="coordinates for each tree")
parser.add_argument("--dlm", type=str, default=".",
                    help="delimeter denoting species")
parser.add_argument("--nodes", action="store_true",
                    help="calculate node heights for a given quartet")
args = parser.parse_args()


def cMono(tree, taxon):
    """Checks if samples are monophyletic in tree
    """
    return(tree.check_monophyly(values=[taxon], target_attr="species")[0])


def getMonophyletic(treelist, quart, winarray):
    """
    """
    p1, p2, p3, p4 = quart
    mtreelist = []
    for i, t in enumerate(treelist):
        if cMono(t, p1) and cMono(t, p2) and cMono(t, p3) and cMono(t, p4):
            Out = t.get_common_ancestor(t.search_nodes(species=quart[-1]))
            t.set_outgroup(Out)
            mtreelist.append(t)
        else:
            winarray[i] = False
#            t2 = t.collapse_lineage_specific_expansions()
#            print(t2)
#            for node in t.split_by_dupes():
#                print(node)
#            ntrees, ndups, sptrees = t.get_speciation_trees()
#            for spt in sptrees:
#                print(spt)
    return(mtreelist, winarray)


def supportFilt(mtreelist, quart, winarray):
    """
    """
    return(mtreelist, winarray)


def loadtrees(treefile, quart, winlist, winarray, nodes, dlm):
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
                treelist.append(t)
    if not winlist:
        winarray = np.ones(len(treelist), dtype=bool)
    mtreelist, winarray = getMonophyletic(treelist, quart, winarray)
    btreelist, winarray = supportFilt(mtreelist, quart, winarray)
    if nodes:
        nh1, nh2 = nodeHeights(btreelist, quart)
    else:
        nh1 = []
        nh2 = []
    return(treelist, winarray, nh1, nh2)


def nodeHeights(mtreelist, quart):
    """Calculates the node heights in a set of trees.

    Parameters
    ------
    trees: ete3 object, returned from function loadtrees
    quart: list, list of topologies to calculate node heights
    windows: list, genome coordinates for which the trees were made

    Returns
    ------
    nh1: float, average value of T1 for topology
    nh2: float, average value of T2 for topology
    """
    print("calculating node heights for quartets: {}".format(quart))
    nh1 = []
    nh2 = []
    for t in mtreelist:
        P1 = t.search_nodes(species=quart[0])
        P2 = t.search_nodes(species=quart[1])
        P3 = t.search_nodes(species=quart[2])
        for p1, p2 in combinations(P1, P2, 2):
            nh2.append(t.get_distance(p1, p2))
        for p1, p3 in combinations(P1, P3, 2):
            nh1.append(t.get_distance(p1, p3))
        for p2, p3 in combinations(P2, P3, 2):
            nh1.append(t.get_distance(p2, p3))
    print("T1: {}, T2: {}".format(np.mean(np.array(nh1)),
                                  np.mean(np.array(nh2))))
    return(nh1, nh2)


def parseWin(windows):
    """
    """
    winlist = []
    with open(windows, 'r') as win:
        win.next()
        for line in win:
            x = line.strip().split()
            winlist.append("{}-{}".format(x[1], x[2]))
    winarray = np.ones(len(winlist), dtype=bool)
    return(winlist, winarray)


def outputTrees(treelist, winlist, winarray, nh1, nh2):
    """
    """
    # mask winlist and treelist using winarray T/F
    return(None)


if __name__ == "__main__":
    quart = args.groups
    if args.windows:
        winlist, winarray = parseWin(args.windows)
    else:
        winlist = []
        winarray = []
    treelist, winarray, nh1, nh2 = loadtrees(args.treefile, quart, winlist,
                                             winarray, args.nodes, args.dlm)
    outputTrees(treelist, winlist, winarray, nh1, nh2)
