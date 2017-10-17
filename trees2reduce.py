#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 12:18:27 2017
@author: stsmall
create pairwise matrix of robison-foulds distance between phylogenetic trees.
"""
import numpy as np
import ete3
import allel
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-t', "--treefile", type=str, required=True,
                    help="treefile in newick, 1 per line")
parser.add_argument('-o', "--outfile", type=str, default="topo.reduced",
                    help="outfile for reduced topologies")
parser.add_argument('-f', "--freq", type=str, default="topo.freq",
                    help="outfile of topology frequency")
parser.add_argument("--topology", action="store_true",
                    help="use only the topology of the tree and not branch"
                    "length information")
parser.add_argument("--pairwise", action="store_true", help="output a large"
                    "pairwise matrix of Robinson-Foulds distances")
args = parser.parse_args()


def removebranchlengths(trees):
    """Remove branch length info from set of newick trees
    """




    return(trees)


def topofreq(trees, out, freq):
    """reduce the symmetric topologies in a set of newick trees using ete3.
    Count the frequency of each topology.
    """



    return(unique_trees)


def rf_distances(unique_trees):
    """Calculate a large pairwise matrix of robinson-foulds distances between
    newick trees.
    """


if __name__ == "__main__":
    trees = args.treefile
    outfile = args.outfile
    freqfile = args.freqfile
    if args.topology:
        "remove branch length info"
        trees = removebranchlengths(trees)
    uniqtrees = topofreq(trees, outfile, freqfile)
    if args.pairwise:
        rf_distances(uniqtrees)
