#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:33:56 2019

@author: scott
"""

from __future__ import print_function
from __future__ import division
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--trees", type=str, required=True, help="trees file")
parser.add_argument("--coords", type=str, help="list of tree files, same order")
parser.add_argument("--clust", type=int, default=100, help="how many loci"
                    " to cluster for astral")
parser.add_argument("--astral_exe", type=str, help="path to astral exec")
parser.add_argument("--scafs", type=str, required=True)
args = parser.parse_args()


def runAstral(treeFile, clust, astralexe):
    """
    """
    tree_list = []
    with open(treeFile) as trees:
        for line in trees:
            tree_list.append(line.strip())
    f = open("astral.tre", 'a')
    start = 0
    step = clust
    end = start + step
    tree_slice = tree_list[start:end]
    while tree_slice:
        t = open("astral_tmp.tre", 'w')
        for t in tree_slice:
            t.write("{}\n".format(t))
        t.close()
        # run astral
        command = "java -jar ~/programs_that_work/ASTRAL/Astral/astral.5.6.1.jar -i astral_tmp.tre -o astral.out -a ../astral.groups"
        proc = subprocess.Popen(command, shell=True)
        with open("astral.out", 'r') as astral:
            for line in astral:
                f.write(line)
        start = end
        end += step
        tree_slice = tree_list[start:end]
    f.close()
    return(None)


def makeWindows(coordList, clust, scaf):
    """
    """
    start_list = []
    end_list = []
    with open(coordList, 'r') as coords:
        for line in coords:
            s, e = line.split("-")
            start_list.append(s)
            end_list.append(e)
    f = open("{}.windows.out".format(scaf), 'w')
    s_ix = 0
    step = clust
    e_ix = s_ix + step
    c = True
    while c:
        try:
            f.write("{}\t{}\{}\n".format(scaf, start_list[s_ix], end_list[e_ix]))
            s_ix = e_ix
            e_ix += step
        except IndexError:
            f.write("{}\t{}\{}\n".format(scaf, start_list[s_ix], end_list[-1]))
            c = False
    f.close()
    return(None)


if __name__ == "__main__":
    runAstral(args.trees, args.clust, args.astral_exe)
    makeWindows(args.coords, args.clust, args.scaf)