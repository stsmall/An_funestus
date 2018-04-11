#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 18:01:25 2018

@author: scott
"""
from __future__ import print_function
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf',
                    help='path to vcf', required=True)
parser.add_argument('-i', '--integer', help='chunk size', default=100000)
args = parser.parse_args()


def chunkVCF(vcfFile, chunk_size):
    """
    """
    sitedict = defaultdict(list)
    header = []
    i = 1
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                header.append(line)
            else:
                j = 1
                while j <= chunk_size:
                    sitedict[i].append(line)
                    try:
                        line = line.next()
                    except EOFError:
                        break
                    j += 1
                i += 1
    for k in sitedict.keys():
        f = open("{}.{}.vcf".format(vcfFile, k), 'w')
        for h in header:
            f.write(h)
        for snp in sitedict[k]:
            f.write(snp)
        f.close()
    return(sitedict)


if __name__ == "__main__":
    chunkVCF(args.vcf, args.integer)
