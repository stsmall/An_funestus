#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 15:07:45 2018
convert LDjump file to recomb map file
@author: stsmall
"""
from __future__ import print_function
from __future__ import division

import argparse
import numpy as np
from math import log

parser = argparse.ArgumentParser()
parser.add_argument("--ld", required=True,
                    help="ldjump file")
parser.add_argument("--Ne", type=int, default=1E6,
                    help="effective size")
parser.add_argument("--mapsize", type=int, required=True,
                    help="mapsize of Chromosome")
args = parser.parse_args()


def readLDjump(LDjumpFile):
    """
    """
    snplist = []
    rholist = []
    with open(LDjumpFile) as rhomap:
        for line in rhomap:
            try:
                x = line.split(",")
                snp = x[0]
                rho = x[1]
                snplist.append(int(snp))
                rholist.append(float(rho))
            except ValueError:
                continue
    return(snplist, rholist)


def makeRecombMap(snplist, rholist, Ne, size):
    """Following the conversion method of Booker et al 2017 in genetics
    """
    poslist = []
    cMMblist = []
    cMlist = []
    cM = 0
    # cumRho = 0
    for i, pos in enumerate(snplist):
        cMMb_avg = 50*log(1/(1-2*(rholist[i]/(4*Ne)))) * 1E6
        poslist.append(pos)
        cMMb_rho = rholist[i] * cMMb_avg  # average rate from Chan 2012
        cMMblist.append(cMMb_rho)  # rate between SNPs
        if i == 0:
            cM += (cMMb_rho * (pos)) / size
            # cumRho += (pos * rholist[i])/size
        else:
            cM += (cMMb_rho * (pos - snplist[i-1])) / size
            # cumRho += (pos-snplist[i-1] * rholist[i])/size
        # print(cumRho)
        cMlist.append(cM)
    return(poslist, cMMblist, cMlist)


if __name__ == "__main__":
    s, r = readLDjump(args.ld)
    poslist, cmMblist, cmlist = makeRecombMap(s, r, args.Ne, args.mapsize)
    with open("cMMb.LD.out", 'w') as cm:
        for i, j in zip(poslist, cmMblist):
            cm.write("{} {}\n".format(i, j))
    with open("cM.LD.out", 'w') as cm:
        for i, j in zip(poslist, cmlist):
            cm.write("{} {}\n".format(i, j))
