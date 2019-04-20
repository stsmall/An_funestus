#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 17:31:59 2019

@author: scott
"""

from bisect import bisect_left as left
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--fdFile", type=str, required=True)
parser.add_argument('-c', "--cMMb", type=str, required=True)
args = parser.parse_args()


def readcMMb(recomb):
    """
    """
    rholist = []
    coordlist = []
    with open(recomb, 'r') as r:
        for line in r:
            coord, rho = line.split()
            rholist.append(float(rho))
            coordlist.append(int(coord))
    return(coordlist, rholist)


def calcAvgRho(fdFile, coordlist, rholist):
    """
    """
    f = open("fdvRho.out", 'w')
    with open(fdFile, 'r') as fd:
        for line in fd:
            chrm, start, end, fd = line.split()
            start_ix = left(coordlist, int(start))
            end_ix = left(coordlist, int(end))
            if end_ix > start_ix:
                rho = np.mean(rholist[start_ix:end_ix])
            if float(fd) < 0:
                fd = 0
            f.write("{}\t{}\t{}\n".format(chrm, rho, fd))
    f.close()
    return(None)


if __name__ == "__main__":
    coordlist, rholist = readcMMb(args.cMMb)
    calcAvgRho(args.fdFile, coordlist, rholist)
            