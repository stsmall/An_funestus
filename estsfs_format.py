#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 16:28:21 2018
Create input files for est-sfs (Keightly and Jackson 2018)

python -i ingroup.counts -o1 outgroup.counts [-o2 outgroup.counts -o3 outgroup.counts]

file format should be from vcftools --counts

@author: stsmall
"""

from __future__ import print_function
from __future__ import division
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--ingroup", type=str, required=True,
                    help="ingroup/focalgroup counts")
parser.add_argument('-n', "--haplotypes", type=int, required=True,
                    help="number of haplotypes in focal group")
parser.add_argument('-ogrpn', "--outgrp_haplotypes", type=int, nargs='+',
                    action="append", required=True,
                    help="number of haplotypes in focal group")
parser.add_argument('-o1', "--outgroup1", type=str, required=True,
                    help="outgroup counts")
parser.add_argument('-o2', "--outgroup2", type=str,
                    help="outgroup counts")
parser.add_argument('-o3', "--outgroup3", type=str,
                    help="outgroup counts")
args = parser.parse_args()


def alleleCounts(fileList):
    """Read in counts print out files
    """
    ancdict = defaultdict(list)
    bporderlist = []
    for f in fileList:
        with open(f, 'r') as ing:
            line = next(ing)
            for line in ing:
                # ACGT
                x = line.split()
                anclist = [0, 0, 0, 0]
                # first allele
                if 'A' in x[4]:
                    anclist[0] += int(x[4].split(":")[1])
                elif 'C' in x[4]:
                    anclist[1] += int(x[4].split(":")[1])
                elif 'G' in x[4]:
                    anclist[2] += int(x[4].split(":")[1])
                elif 'T' in x[4]:
                    anclist[3] += int(x[4].split(":")[1])
                # second allele
                if 'A' in x[5]:
                    anclist[0] += int(x[5].split(":")[1])
                elif 'C' in x[5]:
                    anclist[1] += int(x[5].split(":")[1])
                elif 'G' in x[5]:
                    anclist[2] += int(x[5].split(":")[1])
                elif 'T' in x[5]:
                    anclist[3] += int(x[5].split(":")[1])
                b = "{}_{}".format(x[0], x[1])
                ancdict[b].append(anclist)
                bporderlist.append(b)
    return(ancdict, bporderlist)


if __name__ == "__main__":
    fileList = [args.ingroup, args.outgroup1]
    if args.outgroup2:
        fileList.append(args.outgroup2)
    elif args.outgroup3:
        fileList.append(args.outgroup3)
    haps = args.haplotypes
    haplist = args.outgrp_haplotypes[0]
    ancdict, blist = alleleCounts(fileList)
    t = open("excluded_sites.estsfs.out", 'w')
    with open("est-sfs.data.txt", 'w') as f:
        import ipdb;ipdb.set_trace()
        for k in blist:
            if sum(ancdict[k][0]) != haps:
                t.write("{} {}\n".format(k.split("_")[0], k.split("_")[1]))
            else:
                for i in range(1, len(fileList)):
                    try:
                        ix = ancdict[k][i].index(haplist[i])
                        ancdict[k][i][ix] = 1
                    except ValueError:
                        ancdict[k][i] = [0,0,0,0]
                alleleCt = ["{}".format(",".join(map(str, i))) for i in ancdict[k]]
                f.write("{} {} {}\n".format(k.split("_")[0], k.split("_")[1],
                        " ".join(alleleCt)))
    t.close()
