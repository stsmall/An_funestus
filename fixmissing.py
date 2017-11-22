#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 15:29:18 2017
fixmissing
@author: stsmall
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, action="store",
                    required=True, help="name of infile")
args = parser.parse_args()


def fixmiss(vcffile):
    """
    """
    with open("{}.fixmiss".format(vcffile), 'w') as f:
        with open(vcffile, 'r') as t:
            for line in t:
                if not line.startswith("#"):
                    x = line.strip().split()
                    for i, s in enumerate(x):
                        if ".:." in s:
                            gt = s.split(":")
                            gt[0] = "./."
                            x[i] = ":".join(gt)
                            f.write("{}\n".format("\t".join(x)))
                        else:
                            f.write(line)
    return(None)


if __name__ == "__main__":
    fixmiss(args.vcfFile)
