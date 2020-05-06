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
parser.add_argument('-p1', "--noMigprob", type=float, default=0.10,
                    help="prob of noMig class")
parser.add_argument('-p2', "--migprob", type=float, default=0.90,
                    help="cut off of single class rolling mean")
args = parser.parse_args()


def clusterIntrogressedRegions(InFile, p1, p2):
    """

    """
    clustlist = []
    f = open("{}.bed".format(InFile), 'w')
    with open(InFile, 'r') as filet:
        next(filet)
        for line in filet:
            try:
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
                    if pred == '1':
                        clustlist.append(mig12p)
                        line = next(filet)
                        clustlist.append(float(line.split()[6]))
                        rollmean = np.mean(clustlist)
                        noMigp = float(line.split()[5])
                        while noMigp <= p1 and rollmean >= p2:
                            end = int(line.split()[2])
                            pred += line.split()[4]
                            sites += int(line.split()[3])
                            line = next(filet)
                            clustlist.append(float(line.split()[6]))
                            rollmean = np.mean(clustlist)
                            noMigp = float(line.split()[5])
                        f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, sites, pred, rollmean))
                        clustlist = []
                    elif pred == '2':
                        clustlist.append(mig21p)
                        line = next(filet)
                        clustlist.append(float(line.split()[7]))
                        rollmean = np.mean(clustlist)
                        noMigp = float(line.split()[5])
                        while noMigp <= p1 and rollmean >= p2:
                            end = int(line.split()[2])
                            pred += line.split()[4]
                            sites += int(line.split()[3])
                            line = next(filet)
                            clustlist.append(float(line.split()[7]))
                            rollmean = np.mean(clustlist)
                            noMigp = float(line.split()[5])
                        f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, sites, pred, rollmean))
                        clustlist = []
            except StopIteration:
                if len(clustlist) > 0:
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, sites, pred, rollmean))
                break
    f.close()
    return(None)


if __name__ == "__main__":
    clusterIntrogressedRegions(args.InFile, args.noMigprob, args.migprob)
