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
import matplotlib.pyplot as plt
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


def LoadTrees(treefile, quart, outgroup, dlm):
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


def WriteTrees(treelist):
    """Rewrite reformatted or rerooted trees

    Parameters
    ------
    treelist: list, list of tree objects

    Returns
    ------
    file: file, writes out tree objects in nexus format

    """
    f = open("rooted_trees.nex", 'w')
    for t in treelist:
        f.write("{}\n".format(t.write()))
    f.close()
    return(None)


def cMono(tree, taxon):
    """Checks if samples are monophyletic in tree
    """
    return(tree.check_monophyly(values=taxon, target_attr="species")[0])


def AgeAndSupport(treelist, taxon):
    """Calculates the support and node age if groups in taxon are monophyletic
    """
    taxdict = {}
    for i, tax in enumerate(taxon):
        nodesupport = []
        nodeage = []
        for t in treelist:
            if cMono(t, tax):
                samples = []
                for sp in tax:
                    samples.extend(t.search_nodes(species=sp))
                ancnode = t.get_common_ancestor(samples)
                nodeage.append(ancnode.dist)
                nodesupport.append(ancnode.support)
                # alternate way
#                t2 = t.collapse_lineage_specific_expansions()
#                samples = t2.get_leaf_names()
#                sampletax = [i.split("_")[0] for i in samples]
#                ix = []
#                for species in tax:
#                    ix.append(samples[species.index(sampletax)])
#                nd = t2.get_common_ancestor(ix)
#                nd.dist
#                nd.support
        taxdict[i] = (nodeage, nodesupport)
    return(taxdict)


def AgeStats(taxdict, quart):
    """Plot histogram and stats for branch lengths and support
    """
    # age hist
    for i in taxdict.keys():
        age = taxdict[i][0]
        ma = np.mean(age)  # mean age
        nt = len(age)  # num trees
        print("\n{0}, numtrees {1}, meanage {2:.6f}".format(quart[i], nt, ma))
        label = "".join([j[0] for j in quart[i]])
        plt.hist(age, bins=30, alpha=0.5, label=label)
    plt.legend(loc='upper right')
    plt.ylabel("node age")
    plt.savefig("NodeHeight.pdf")
    plt.clf()
    return(None)


def SupportStats(taxdict, quart):
    # support hist
    for i in taxdict.keys():
        supp = taxdict[i][1]
        nt = len(supp)  # num trees
        ms = np.mean(supp)  # mean support
        print("\n{0}, numtrees {1}, meansupport {2:.2f}".format(quart[i], nt, ms))
        label = "".join([j[0] for j in quart[i]])
        plt.hist(supp, bins=30, alpha=0.5, label=label)
    plt.legend(loc='upper right')
    plt.ylabel("support")
    plt.savefig("NodeSupport.pdf")
    plt.clf()
    return(None)


def SupportFilt(treelist, quart):
    """Collapses nodes or make polytomy with support below threshold
    """
    return(None)


def ParseWin(windows):
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


def TreesTable(treelist, winlist, winarray, nh1, nh2):
    """Print windowed statistics for trees, also print corresponding trees
    """
    # remove trees with low support
    # write rerooted trees
    # mask winlist and treelist using winarray T/F
    return(None)


if __name__ == "__main__":
    quart = args.groups
    treelist = LoadTrees(args.treefile, quart, args.outgroup, args.dlm)
    WriteTrees(treelist)
    taxdict = AgeAndSupport(treelist, quart)
    AgeStats(taxdict, quart)
    SupportStats(taxdict, quart)
