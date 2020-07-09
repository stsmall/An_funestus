#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 13:23:11 2018
pruneTips.py -t TREE.tre -g groups -n 1 -dlm
@author: scott
"""
import numpy as np
from ete3 import PhyloTree
import re
from tqdm import tqdm, trange
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', "--treefile", type=str, required=True,
                    help="treefile in newick, 1 per line")
parser.add_argument('-g', "--groups", nargs='+', action='append',
                    required=True, help="quartet of species to calculate,"
                    " assumes form: P1 P2 P3. can be given multiple times")
parser.add_argument("--dlm", type=str, default="_",
                    help="delimeter denoting species")
parser.add_argument("--rand", action="store_true", help="choose individual at"
                    "random note you will have to format leaf names")
parser.add_argument("--rename", action="store_true", help="rename tips")
args = parser.parse_args()


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return(i + 1)


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
    pbar = tqdm(total=file_len(treeFile))
    with open(treeFile, 'r') as newick:
        for line in newick:
            pbar.update(1)
            if not line.startswith("NA"):
                t = PhyloTree(line)
                t.set_species_naming_function(lambda node: node.name.split(dlm)[0])
                treelist.append(t)
    pbar.close()
    return(treelist)


def pruneTips(treelist, species, rand, topo=True, ntaxa=1):
    """Prune trees so only n taxa remain from each of species
    """
    print(f"pruning tips ...")
    pbar = tqdm(total=len(treelist))
    if rand:  # keep 1 or more
        splist = []
        for tax in species:
            splist.append(treelist[0].search_nodes(species=tax))
        for t in treelist:
            pbar.update(1)
            nodes = [np.random.choice(g, ntaxa) for g in splist]
            nds = np.concatenate(nodes).ravel()
            t.prune(nds, preserve_branch_length=topo)
    else:  # if only give a single individual
        for t in treelist:
            pbar.update(1)
            t.prune(species, preserve_branch_length=topo)
    pbar.close()
    return(treelist)


def WriteTrees(treelist, rand, rename):
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
        if rename:
            t2 = re.sub(r'(\.[A-Z])\w+', '', t.write())
            f.write("{}\n".format(t2))
        else:
            f.write("{}\n".format(t.write()))
    f.close()
    return(None)


if __name__ == "__main__":
    species = args.groups
    treelist = LoadTrees(args.treefile, args.dlm)
    pruneTips(treelist, species[0], args.rand)
    WriteTrees(treelist, args.rand, args.rename)
