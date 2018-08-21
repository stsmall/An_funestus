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

parser = argparse.ArgumentParser()
parser.add_argument("--ld", required=True,
                    help="ldjump file")
args = parser.parse_args()


def readLDjump(LDjumpFile):
    """
    """
    snplist = []
    rholist = []
    with open(LDjumpFile) as rhomap:
        for line in rhomap:
            x = line.split()
            snp = x[0]
            rho = x[1]
            snplist.append(int(snp))
            rholist.append(float(rho))
    return(snplist, rholist)


def makeRecombMap(snplist, rholist):
    """Following the conversion method of Booker et al 2017 in genetics
    """
    poslist = []
    cmlist = []
    crho = 0
    for i, pos in enumerate(snplist):
        try:
            crho += rholist[i] * (pos - snplist[i-1])
        except IndexError:
            crho += rholist[i] * (pos - 0)
        poslist.append(pos)
        cmlist.append(crho/pos)
    return(poslist, cmlist)


if __name__ == "__main__":
    s, r = readLDjump(args.ld)
    poslist, cmlist = makeRecombMap(s, r)
    with open("cM.LD.out", 'w') as cm:
        for i, j in zip(poslist,cmlist):
            cm.write("{} {}\n".format(i, j))
