#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 13:26:50 2018

@author: scott
"""

from __future__ import print_function
import argparse
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--transferFile", action='store', required=True,
                    help="Chromosome and position transfer table")
args = parser.parse_args()


def choose_duplicates(transferFile):
    """Liftover file seems to contain duplicate mappings
    """
    paralog = open("{}.paralogs".format(transferFile), 'w')
    transdict = defaultdict(dict)
    with open(transferFile, 'r') as transfer:
        for line in transfer:
            x = line.strip().split()
            chrom = x[4]
            orient_r = x[1]  # this should be positive for the ref
            if orient_r != "+":
                raise Exception("Reference/Query Orientation should be + orientation")
            pos_e = x[7]
            try:
                transdict[chrom][pos_e]
                # duplicate
                first_line = transdict[chrom][pos_e]
                if chrom in first_line.split()[0]:
                    pass
                elif chrom in line.split()[0]:
                    transdict[chrom][pos_e] = line
                else:
                    pass
                paralog.write("{}\n{}".format(transdict[chrom][pos_e], line))
           except KeyError:
                transdict[chrom][pos_e] = line
    paralog.close()
    return(transdict)


if __name__ == "__main__":
    transdict = choose_duplicates(args.transferFile)
    f = open("{}.nodups".format(args.transferFile), 'w')
    for k in transdict.keys():
        for p in transdict[k].keys():
            f.write("{}".format(transdict[k][p]))
    f.close()
