#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 15:44:31 2018

@author: scott
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--InFile", required=True, help="filet output file")
parser.add_argument('-p1', "--prob1", type=float, default=0.05,
                    help="reject no-migration class")
parser.add_argument('-p2', "--prob2", type=float, default=0.90,
                    help="accept migration class")
parser.add_argument("--slide", action="store_true")
args = parser.parse_args()


def clusterIntrogressedRegions(InFile, p1, p2, slide):
    """Split introgressed regions to file."""
    c1 = open(f"{InFile}.m12.bed", 'w')
    c2 = open(f"{InFile}.m21.bed", 'w')
    c3 = open(f"{InFile}.m.bed", 'w')
    c4 = open(f"{InFile}.no.bed", 'w')
    with open(InFile, 'r') as filet:
        next(filet)
        for line in filet:
            x = line.split()
            chrom = x[0]
            start = int(x[1])
            end = int(x[2])
            if slide:
                sites = x[-1]
            else:
                sites = int(x[3])
            pred = x[4]
            noMigp = float(x[5])
            mig12p = float(x[6])
            mig21p = float(x[7])
            if noMigp <= p1:
                if mig12p >= p2:
                    c1.write(f"{chrom}\t{start}\t{end}\t{sites}\t{pred}\t{mig12p}\n")
                elif mig21p >= p2:
                    c2.write(f"{chrom}\t{start}\t{end}\t{sites}\t{pred}\t{mig21p}\n")
                else:
                    # introgression but not high confidence for direction
                    c3.write(f"{chrom}\t{start}\t{end}\t{sites}\t{pred}\t{noMigp}\t{mig12p}\t{mig21p}\n")
            else:
                c4.write(f"{chrom}\t{start}\t{end}\t{sites}\t{pred}\t{noMigp}\n")
    c1.close()
    c2.close()
    c3.close()
    c4.close()
    return(None)


if __name__ == "__main__":
    clusterIntrogressedRegions(args.InFile, args.prob1, args.prob2. args.slide)
