#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 17:25:00 2018

@author: scott
"""

from __future__ import print_function
from __future__ import division
import sys
import bisect
import argparse
# check versions
assert sys.version_info[0] >= 3

parser = argparse.ArgumentParser()
parser.add_argument('-vcf', "--vcfFile", required=True,
                    help="vcf for thinning")
parser.add_argument('-gff', "--gffFile", required=True,
                    help="gff with coords")
parser.add_argument('-t', "--thin", type=int, default=5000,
                    help="distance from gene")
args = parser.parse_args()


def readgff(gffFile):
    """
    """
    gffdict = {}
    with open(gffFile, 'r') as gff:
        for line in gff:
            if line.startswith("#"):
                pass
            else:
                x = line.split()
                if x[2] == 'gene':
                    chrom = x[0]
                    try:
                        gffdict[chrom][0].append(int(x[3]))
                        gffdict[chrom][1].append(int(x[4]))
                    except KeyError:
                        gffdict[chrom] = ([], [])
                        gffdict[chrom][0].append(int(x[3]))
                        gffdict[chrom][1].append(int(x[4]))
    return(gffdict)


def thinVcf(gffdict, vcfFile, thin):
    """
    """
    f = open("vcf.out", 'w')
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                pos = int(x[1])
                chrom = x[0]
                import ipdb;ipdb.set_trace()
                s = bisect.bisect_right(pos, gffdict[chrom][0])
                e = bisect.bisect_right(pos, gffdict[chrom][1])
                if s-1 == e:
                    # in a gene
                    pass
                else:
                    if pos + thin < gffdict[chrom][0][s] and pos - thin > gffdict[chrom][1][s-1]:
                        f.write(line)
    return(None)


if __name__ == "__main__":
    gffdict = readgff(args.gffFile)
    thinVcf(gffdict, args.vcfFile, args.thin)
