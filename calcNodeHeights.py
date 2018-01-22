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
    agelist = []
    supportlist = []
    for tax in taxon:
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
            else:
                nodeage.append(np.nan)
                nodesupport.append(np.nan)
        agelist.append(nodeage)
        supportlist.append(nodesupport)
    taxdict["age"] = agelist
    taxdict["support"] = supportlist
    return(taxdict)


def AgeStats(taxdict, quart):
    """Plot histogram and stats for branch lengths and support
    """
    # age hist
    for pop, age in enumerate(taxdict['age']):
        ma = np.nanmean(age)  # mean age
        nt = np.count_nonzero(~np.isnan(age))  # num trees
        print("\n{0}, numtrees {1}, meanage {2:.6f}".format(quart[pop], nt, ma))
        label = "".join([j[0] for j in quart[pop]])
        age = np.array(age)
        age = age[~np.isnan(age)]
        plt.hist(age, bins=30, alpha=0.5, label=label)
    plt.legend(loc='upper right')
    plt.ylabel("node age")
    plt.savefig("NodeHeight.pdf")
    plt.clf()
    return(None)


def SupportStats(taxdict, quart):
    # support hist
    for pop, supp in enumerate(taxdict['support']):
        ms = np.nanmean(supp)  # mean age
        nt = np.count_nonzero(~np.isnan(supp))  # num trees
        print("\n{0}, numtrees {1}, meansupport {2:.2f}".format(quart[pop], nt, ms))
        label = "".join([j[0] for j in quart[pop]])
        supp = np.array(supp)
        supp = supp[~np.isnan(supp)]
        plt.hist(supp, bins=30, alpha=0.5, label=label)
    plt.legend(loc='upper right')
    plt.ylabel("support")
    plt.savefig("NodeSupport.pdf")
    plt.clf()
    return(None)


def WindowStats(windows, taxdict, quart):
    """
    """
    winlist = []
    f = open("window_stats.tsv", 'w')
    with open(windows, 'r') as win:
        header = win.next()
        for line in win:
            winlist.append(line.rstrip("\n"))
    label = []
    for tax in quart:
        label.append("".join([j[0] for j in tax]))
    f.write("Summary\t{}\t{}\n".format(header.rstrip("\n"), "\t".join(label)))
    nage = zip(*taxdict['age'])
    nsupp = zip(*taxdict['support'])
    try:
        for i, age in enumerate(nage):
            f.write("NodeHeight\t{}\t{}\n".format(winlist[i], "\t".join(map(str, age))))
            f.write("NodeSupport\t{}\t{}\n".format(winlist[i], "\t".join(map(str, nsupp[i]))))
    except TypeError:
        import ipdb;ipdb.set_trace()
    f.close()
    return(None)


if __name__ == "__main__":
    quart = args.groups
    treelist = LoadTrees(args.treefile, quart, args.outgroup, args.dlm)
    WriteTrees(treelist)
    taxdict = AgeAndSupport(treelist, quart)
    AgeStats(taxdict, quart)
    SupportStats(taxdict, quart)
    if args.windows:
        WindowStats(args.windows, taxdict, quart)

