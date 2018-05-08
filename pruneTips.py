#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 13:23:11 2018
pruneTips.py -t TREE.tre -g groups -n 1 -dlm
@author: scott
"""

from __future__ import print_function
from __future__ import division
import numpy as np
from ete3 import PhyloTree
import re

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', "--treefile", type=str, required=True,
                    help="treefile in newick, 1 per line")
parser.add_argument('-g', "--groups", nargs='+', action='append',
                    required=True, help="quartet of species to calculate,"
                    " assumes form: P1 P2 P3. can be given multiple times")
parser.add_argument('-n', "--numberTips", type=int, default=1,
                    help="number of tips to have after pruning")
parser.add_argument("--dlm", type=str, default="_",
                    help="delimeter denoting species")
parser.add_argument("--rand", action="store_true", help="choose individual at"
                     "random note you will have to format leaf names")
args = parser.parse_args()


def LoadTrees(treeFile, dlm):
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
    with open(treeFile, 'r') as newick:
        for line in newick:
            if not line.startswith("NA"):
                t = PhyloTree(line)
                t.set_species_naming_function(lambda node: node.name.split(dlm)[0])
                treelist.append(t)
    return(treelist)


def pruneTips(treelist, species, n, rand, topo=False):
    """Prune trees so only n taxa remain from each of species
    """
    if rand:
        splist = []
        for tax in species:
            splist.append(treelist[0].search_nodes(species=tax))
        for t in treelist:
            while error:
                try:
                    nodes = [np.random.choice(g, n) for g in splist]
                    nds = np.concatenate(nodes).ravel()
                    t.prune(nds, preserve_branch_length=topo)
                    error = False
                except:
                    error = True
    else:
        for t in treelist:
            t.prune(species, preserve_branch_length=topo)
    return(treelist)


def WriteTrees(treelist):
    """Rewrite reformatted or rerooted trees

    Parameters
    ------
    treelist: list, list of tree objects

    Returns
    ------
    file: file, writes out tree objects in nexus format

    """
    f = open("trees.rp.nex", 'w')
    for t in treelist:
        t2 = re.sub(r'_([0-9]|[A-Z])\w+', '', t.write())
        f.write("{}\n".format(t2))
    f.close()
    return(None)


if __name__ == "__main__":
    species = args.groups
    treelist = LoadTrees(args.treefile, args.dlm)
    pruneTips(treelist, species[0], args.numberTips, args.rand)
    WriteTrees(treelist)
