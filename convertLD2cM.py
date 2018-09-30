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
from math import log

parser = argparse.ArgumentParser()
parser.add_argument("--ld", required=True, help="ldjump file")
parser.add_argument("--chrom", type=str, required=True, help="chrom map size")
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
        cMlist.append(cM)
    return(poslist, cMMblist, cMlist)


def makeRecombMapBooker(snplist, rholist, chrom):
    """Following the conversion method of Booker et al 2017 in genetics
    """
    if chrom == "X":
        mapsize = 44.7
    elif chrom == "3R":
        mapsize = 90.9
    elif chrom == "3L":
        mapsize = 89.1
    elif chrom == "2L":
        mapsize = 63.2
    elif chrom == "2R":
        mapsize = 94.8
    poslist = []
    rhocum = []
    cMlist = []
    cMMblist = []
    for i, pos in enumerate(snplist):
        if i == 0:
            rhoTemp = (rholist[i] * (pos))
        else:
            rhoTemp = (rholist[i] * (pos - snplist[i-1]))
        if i == 0:
            rhocum.append(rhoTemp)
        else:
            rhocum.append(rhocum[-1] + rhoTemp)
        poslist.append(pos)
    for i, j in enumerate(rhocum):
        cMperSNP = (j / rhocum[-1])
        cMlist.append(cMperSNP)
        cMMblist.append(((cMlist[i] - cMlist[i-1])*mapsize) / ((snplist[i] - snplist[i-1])/1E6))
    return(poslist, cMMblist, cMlist)


if __name__ == "__main__":
    s, r = readLDjump(args.ld)
#    poslist, cmMblist, cmlist = makeRecombMap(s, r, args.Ne, args.mapsize)
    poslist, cmMblist, cmlist = makeRecombMapBooker(s, r, args.chrom)
    with open("cMMb.LD.out", 'w') as cm:
        for i, j in zip(poslist, cmMblist):
            cm.write("{} {}\n".format(i, j))
    with open("cM.LD.out", 'w') as cm:
        for i, j in zip(poslist, cmlist):
            cm.write("{} {}\n".format(i, j))
