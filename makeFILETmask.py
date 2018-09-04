#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 15:44:31 2018

@author: scott
"""
from __future__ import division
from __future__ import print_function
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--InFile", required=True, help="fasta input file")
parser.add_argument('-w', "--window_size", type=int, default=10000,
                    help="window size for calculation")
parser.add_argument('-n', "--nLength", type=int, default=20, help="length of "
                    "Ns to cosider for masking")
args = parser.parse_args()


def buildMaskFile(InFile, window, nLength):
    """
    """
    r = re.compile(r'N{{{0},}}'.format(nLength))
    f = open("{}.mask".format(InFile), 'w')
    with open(InFile, 'r') as fasta:
        for line in fasta:
            if line.startswith(">"):
                line = next(fasta)
                seq = line.strip()
                start = 0
                end = start + window
                step = window
                while end < len(seq):
                    seqR = seq[start:end]
                    cord = [(m.start(), m.end()) for m in re.finditer(r, seqR)]
                    for i, j in cord:
                        f.write("0 {} {}\n\n//\n\n".format(i/window, j/window))
                    start += step
                    end += step
                break
    f.close()
    return(None)


if __name__ == "__main__":
    buildMaskFile(args.InFile, args.window_size, args.nLength)
