#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 12:48:37 2018

@author: scott
"""

from __future__ import print_function
from __future__ import division
import numpy as np
from itertools import combinations
import argparse
import ipdb

parser = argparse.ArgumentParser()
parser.add_argument('-d', "--distFile", type=str, required=True,
                    help="treefile in newick, 1 per line")
parser.add_argument('-w', "--windows", type=str, required=True,
                    help="coordinates for each tree")
parser.add_argument("--dlm", type=str, default="_",
                    help="delimeter denoting species")
parser.add_argument("-o", "--outgroup", type=str,
                    help="outgroup species for rooting")
args = parser.parse_args()


def getFirstSet(distFile):
    """
    """
    with open(distFile, 'r') as d:
        for line in d:
            if line.strip().isdigit():
                digit = line.strip()
                line = next(d)
                vert = []
                while line.rstrip("\n") != '':
                    x = line.split()
                    vert.append(x[0].split("_")[0])
                    line = next(d)
                if line.rstrip("\n") == '':
                    return(vert, digit)


def countBlocks(distFile):
    """
    """
    count = 0
    with open(distFile, 'r') as d:
        for line in d:
            if line.strip().isdigit():
                count += 1
    return(count)


def parseDistmat(distFile):
    """
    """
    count = countBlocks(distFile)
    vert, digit = getFirstSet(distFile)
    pwmatrix = np.zeros([count, int(digit), int(digit)])
    i = 0
    try:
        with open(distFile, 'r') as d:
            for line in d:
                if line.strip().startswith(digit):
                    line = d.next()
                    x = []
                    while line.rstrip("\n") != '':
                        x.append(line.split()[1:])
                        line = d.next()
                    pwmatrix[i] = x
                    i += 1
    except StopIteration:
        pwmatrix[i] = x
        return(pwmatrix, (vert, digit, count))


def parseWindows(windowFile):
    """
    """
    winlist = []
    with open(windowFile, 'r') as win:
        header = win.next()
        for line in win:
            winlist.append(line.rstrip("\n"))
    return(header, winlist)


def printDistmat(pwmatrix, vdc, dlm, outgroup, header, winlist):
    """
    """
    f = open('distwindows.tsv', 'w')
    vert, digit, count = vdc
    spcomb = []
    for p in combinations(vert, 2):
        spcomb.append("{}_{}".format(p[0], p[1]))
    f.write("{}\t{}\n".format(header.rstrip("\n"), "\t".join(spcomb)))
    ipdb.set_trace()



    f.close()
    return(None)


if __name__ == "__main__":
    header, winlist = parseWindows(args.windows)
    pwmatrix, vdc = parseDistmat(args.distFile)
    printDistmat(pwmatrix, vdc, args.dlm, args.outgroup, header, winlist)
