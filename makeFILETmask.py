#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 15:44:31 2018
usage: python makeFILETmask.py -f FILE -n 100 -w 10000
requires that fasta is in 1 line and not at character breaks
>2L
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
parser.add_argument('-n', "--nLength", type=int, default=100, help="length of "
                    "Ns to cosider for masking")
parser.add_argument('-m', "--skipMask", type=float, default=0.50)
parser.add_arguement('-n', "--numberMask", type=int, help="how many mask"
                     "lines")
args = parser.parse_args()


def buildMaskFile(InFile, window, nLength, skipMask, numb):
    """
    """
    n = 0
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
                    if (seqR.count('N') / window) <= 0.50:
                        cord = [(m.start(), m.end()) for m in re.finditer(r, seqR)]
                        if cord:
                            for i, j in cord:
                                f.write("0 {} {}\n".format(i/window, j/window))
                        else:
                            f.write("0 0 0\n")
                        f.write("\n//\n\n")
                    start += step
                    end += step
                    n += 1
                    if n == numb:
                        break
                break
    f.close()
    return(None)


if __name__ == "__main__":
    buildMaskFile(args.InFile,
                  args.window_size,
                  args.nLength,
                  args.skipMask,
                  args.numberMask)
