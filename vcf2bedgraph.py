#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 14:41:17 2017
ref2bedgraphy.py vcf
@author: stsmall
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf file')
args = parser.parse_args()


def ref2bg(vcfFile):
    """
    """
    f = open("{}.bedgraph".format(vcfFile), 'w')
    f.write("track type=bedGraph\n")
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if not line.startswith("#"):
                x = line.strip().split()
                chrom = x[0]
                pos_e = x[1]
                pos_s = int(pos_e) - 1
                ref_a = x[3]
                f.write("{}\t{}\t{}\t{}\n".format(chrom, pos_s, str(pos_e),
                                                  ref_a))
    f.close()
    return(None)


if __name__ == '__main__':
    ref2bg(args.INvcf)
