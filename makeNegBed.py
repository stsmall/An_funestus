#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 14:24:01 2018
python makeNegBed.py 2L.gff 2L
alternatively: bedtools flank +-2000, bedtools complement, bedtools windows
@author: scott
"""

import sys

distance = 2000
length = 1000
cdsdict = {}
i = 0
f = open("{}.neg.bed".format(sys.argv[1]), 'w')
with open(sys.argv[1], 'r') as gff:
    for line in gff:
        if not line.startswith("#"):
            x = line.split()
            chrom = x[0]
            start = x[3]
            end = x[4]
            if "gene" in x[2]:
                cdsdict["cds_" + str(i)] = [int(start)-1, int(end)]
                i += 1

for i, key in enumerate(cdsdict.keys()):
    k = "cds_" + str(i)
    if i == 0:
        end = 0
        next_start = cdsdict[k][0]
    else:
        try:
            end = cdsdict[k][1]
            next_start = cdsdict["cds_" + str(i+1)][0]
        except KeyError:
            break
    Sstart = end + distance
    Send = Sstart + length
    if next_start - Send > distance:
        f.write("{}\t{}\t{}\n".format(sys.argv[2], int(Sstart)+1, int(Send)))
f.close()
