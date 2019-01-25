# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 12:18:27 2017
@author: stsmall
create pairwise matrix of robison-foulds distance between phylogenetic trees.
"""
from __future__ import print_function
from __future__ import division

import numpy as np
from ete3 import Tree
import re
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-t', "--treefile", type=str, required=True,
                    help="treefile in newick, 1 per line")
parser.add_argument('-w', "--windows", type=str,
                    help="coordinates for each tree")
parser.add_argument('-r', "--reftree", type=str,
                    help="RF distance from a fixed topology")
parser.add_argument('-o', "--outgroup", type=str,
                    help="outgroup seq")
parser.add_argument("--threshold", type=int, default=0,
                    help="trees less than this value are considered the same")
parser.add_argument("--topology", action="store_true",
                    help="use only the topology of the tree and not branch"
                    "length information")
parser.add_argument("--pairwise", action="store_true", help="output a large"
                    "pairwise matrix of Robinson-Foulds distances")
parser.add_argument('-g', "--groups", nargs='+', help="list of groups to "
                     "check for monophyly in the trees")
parser.add_argument("--nodeheights", nargs='+', help="calculate node heights"
                     "for a quartet")
args = parser.parse_args()

# TODO: check monophyly for list of groups
# TODO: calculate node heights for specific topologies


def removebl(newick_tree):
    """Remove branch length info from set of newick trees. This simplifies the
    comparisons when searching for non-monophyletic trees

    Parameters
    ------
    newick_tree: string, single tree in newick format

    Returns
    -----
    newick_tree: string, single tree in newick format with no branch length
        info
    """
    print("removing branch lengths ... ")
    # use regex to replace branch lengths by empty
    topolist2 = [re.sub(r'(:\-?[0-9]e?\.?\-?[0-9]*)', '', t) for t in treelist]
    topolist = [re.sub(r'e-[0-9]*', '', t) for t in topolist2]
    return(topolist)


def loadtrees(newickfile, topo, outgroup):
    """Reads and stores phylogenetic trees from a file

    Parameters
    ------
    newickfile: file, file of newick trees, 1 per line
    topo: bool, use only the topology and not branch information
    outgroup: string, outgroup species

    Returns
    ------
    treelist: obj, ete3 object of trees

    """
    print("loading trees...")
    treelist = []
    if topo:
        topolist = []
        with open(newickfile, 'r') as t:
            for line in t:
                topolist.append(line.strip())
        topolist = removebl(topolist)
        for tree in topolist:
            t1 = Tree(tree)
            if outgroup:
                t1.set_outgroup(outgroup)
            treelist.append(t1)
    else:
        with open(newickfile, 'r') as t:
            for line in t:
                if not line.startswith("NA"):
                    t1 = Tree(line)
                    if outgroup:
                        t1.set_outgroup(outgroup)
                    treelist.append(t1)
    return(treelist)


def topofreq(trees, rfthresh):
    """reduce the symmetric topologies in a set of newick trees using ete3.
    Count the frequency of each topology. Return a set of unique topologies

    Parameters
    ------
    trees: obj, ete3 obj from loadtrees
    rfthresh: float, threshold where trees with a R-F distance less than this
        are considered identical and collapsed

    Returns
    ------
    utrees: list, list of unique trees

    """
    print("reducing trees...")
    uniqfreq = {}
    treearray = np.array(treelist)
    while len(treearray) > 1:
        tree1, treearray = treearray[-1], treearray[:-1]
        rflist = []
        for tree2 in treearray:
            try:
                rf = tree1.robinson_foulds(tree2)
            except Exception as e:
                if "unrooted" in str(e):
                    rf = tree1.robinson_foulds(tree2, unrooted_trees=True)
                    print("setting as unrooted, consider using --outgroup")
                else:
                    import ipdb;ipdb.set_trace()
            rflist.append(rf[0])
        rfarray = np.array(rflist)
        uniqfreq[tree1] = np.count_nonzero(rfarray <= rfthresh) + 1
        # shorten the list
        symtrees = rfarray > rfthresh
        treearray = treearray[symtrees]
        print(len(treearray))
    print("writing output files...")
    # keep lists ordered
    uniqtrees = []
    freqlist = []
    for kt in uniqfreq:
        uniqtrees.append(kt)
        freqlist.append(uniqfreq[kt])
    # write out topologies
    topofile = open("trees.topo", 'w')
    for t in uniqtrees:
        topofile.write("{}\n".format(t.write(format=9)))
    topofile.close()
    # write out frequencies
    with open("trees.freq", 'w') as freq:
        for f in freqlist:
            freq.write("{}\n".format(f))
    return(uniqtrees)


def pwdistance(utrees):
    """Calculate a large pairwise matrix of robinson-foulds distances between
    newick trees.

    Parameters
    ------
    utrees: list, topologies stored as a string

    """
    print("calculating pairwise distances")
    pwmat = np.zeros([len(uniqtrees), len(uniqtrees)])
    for i, x in enumerate(uniqtrees):
        for j, y in enumerate(uniqtrees):
            pwmat[i, j] = x.robinson_foulds(y)[0]
    np.savetxt("rf.pwmatrix.csv", pwmat, delimiter=',', fmt='%1.2f')
    return(None)


def refdistance(trees, reftree, coords, outgroup):
    """Calculates the RF distance between a reference topology and set of trees
    in windows along a chromosome

    Parameters
    ------
    trees: obj, obj from loadtrees
    reftree: string, reference topology to check distance against
    coords: list, list of genome coordinates for each tree
    outgroup: string, name of outgroup

    Returns
    ------
    file

    """
    print("getting RF distances from ref")
    with open(reftree, 'r') as tt:
        for line in tt:
            reftree = Tree(line)
    if outgroup:
        reftree.set_outgroup(outgroup)
    rfnorm = []
    for t in trees:
        try:
            rf = t.compare(reftree)
        except Exception as e:
            if "unrooted" in str(e):
                rf = t.compare(reftree, unrooted_trees=True)
                print("setting as unrooted, consider using --outgroup")
            else:
                import ipdb;ipdb.set_trace()
        rfnorm.append(rf["norm_rf"])
    with open("rf.reftree", 'w') as f:
        if coords:
            f.write("start\tstop\trfnorm\n")
            reflist = zip(coords, rfnorm)
            for x, y in reflist:
                f.write("{}\t{}\t{}\n".format(x.split("-")[0], x.split("-")[1], y))
        else:
            for y in rfnorm:
                f.write("{}\n".format(y))
    return(None)


if __name__ == "__main__":
    trees = args.treefile
    windows = args.windows
    coords = []
    if windows:
        with open(args.windows, 'r') as w:
            for line in w:
                if not line.startswith("#"):
                    x = line.strip().split()
    #                s = re.search(r'[S-s]'x.index("Start")
    #                e = x.index("End")
                    coord = "{}-{}".format(x[1], x[2])
                    coords.append(coord)
    treelist = loadtrees(trees, args.topology, args.outgroup)
    uniqtrees = topofreq(treelist, args.threshold)
    if args.pairwise:
        pwmat = pwdistance(uniqtrees)
    if args.reftree:
        refdistance(treelist, args.reftree, coords, args.outgroup)
