#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 15:17:54 2017
fillrefinfb.py
@author: stsmall
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf file')
args = parser.parse_args()


def fixref(vcfFile):
    """
    """
    f = open("{}.reffill".format(vcfFile), 'w')
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if line.startswith("##"):
                f.write(line)
            elif line.startswith("#CHROM"):
                f.write(line)
            else:
                x = line.strip().split()
                if "./." not in x[9].split(":")[0]:
                    for i, s in enumerate(x[10:]):
                        gt = s.split(":")
                        if "./." in gt[0]:
                            gt[0] = "0/0"
                            x[i + 10] = ":".join(gt)
                    f.write("{}\n".format("\t".join(x)))
                else:
                    f.write(line)
    f.close()
    return(None)


if __name__ == '__main__':
    fixref(args.INvcf)
