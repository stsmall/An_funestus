#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 13:42:37 2018

@author: scott
"""
from __future__ import print_function
from __future__ import division
import argparse
from ete3 import Tree
parser = argparse.ArgumentParser()
parser.add_argument('-t', "--treefile", type=str, required=True,
                    help="treefile in newick, 1 per line")
parser.add_argument('-o', "--outgroup", type=str)
args = parser.parse_args()


def makeUltra(treeFile, outgroup):
    """Make a tree ultrametric
    """
    print("loading trees...")
    treelist = []
    with open(treeFile, 'r') as newick:
        for line in newick:
            if not line.startswith("NA"):
                t = Tree(line)
                if outgroup:
                    t.set_outgroup("outgroup")
                t.convert_to_ultrametric()
                treelist.append(t)
    return(treelist)


def writeTrees(treelist):
    """Rewrite reformatted or rerooted trees

    Parameters
    ------
    treelist: list, list of tree objects

    Returns
    ------
    file: file, writes out tree objects in nexus format

    """
    f = open("ultra_trees.nex", 'w')
    for t in treelist:
        f.write("{}\n".format(t.write()))
    f.close()
    return(None)


if __name__ == "__main__":
    tl = makeUltra(args.treefile, args.outgroup)
    writeTrees(tl)
#    t=Tree("(rivulorum:0.150148,(((funestuscf:0.007258,funestus:0.001803)"
#          "1:0.000495,vaneedeni:0.003145)1:0.000816,(longipalpusC:0.001316"
#           ",parensis:0.004594)1:0.008818)1:0.150148);")
