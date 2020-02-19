#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 15:44:31 2018

@author: scott
"""

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--InFile", required=True, help="filet output file")
parser.add_argument('-p1', "--prob1", type=float, default=0.10,
                    help="prob of noMig class")
parser.add_argument('-p2', "--prob2", type=float, default=0.90,
                    help="cut off of single class rolling mean")
args = parser.parse_args()


def clusterIntrogressedRegions(InFile, p1, p2):
    """

    """
    c1 = open("{}.m12.bed".format(InFile), 'w')
    c2 = open("{}.m21.bed".format(InFile), 'w')
    c3 = open("{}.m.bed".format(InFile), 'w')
    with open(InFile, 'r') as filet:
        next(filet)
        for line in filet:
            x = line.split()
            chrom = x[0]
            start = int(x[1])
            end = int(x[2])
            sites = int(x[3])
            pred = x[4]
            noMigp = float(x[5])
            mig12p = float(x[6])
            mig21p = float(x[7])
            if noMigp <= p1:
                if pred == '1' and mig12p >= p2:
                    c1.write(f"{chrom}\t{start}\t{end}\t{sites}\t{pred}\n")
                elif pred == '2' and mig21p >= p2:
                    c2.write(f"{chrom}\t{start}\t{end}\t{sites}\t{pred}\n")
                else:
                    c3.write(f"{chrom}\t{start}\t{end}\t{sites}\t{pred}\t{noMigp}\t{mig12p}\t{mig21p}\n")
    c1.close()
    c2.close()
    c3.close()
    return(None)

if __name__ == "__main__":
    clusterIntrogressedRegions(args.InFile, args.prob1, args.prob2)
