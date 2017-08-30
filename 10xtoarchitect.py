#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 18:33:26 2016
creates input for bam2containmnet.py
consumes a sam file mapped with bwa and option -X and a barcode file a list
of unique barcodes in the mapped sam
@author: scott
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('sam', metavar="sam", type=str,
                    help='path to mapped sam file')
parser.add_argument('-b', '--barcode', type=str, required=True,
                    help="barcode file")
args = parser.parse_args()


def makebamwellID(sam, barcode):
    """
    """
    bc = {}
    wellid = 1
    with open(barcode, 'r') as code:
        for line in code:
            bc[line.strip()] = wellid
            wellid += 1
    f = open(sam+".bc", 'w')
    with open(sam, 'r') as samfile:
        for line in samfile:
            if line.startswith("@"):
                f.write(line)
            else:
                x = line.split()
                if "BX" in x[-1]:
                    ids = x[-1].split(":")[-1]
                    x[0] = "well{}_".format(bc[ids]) + x[0]
                    f.write("\t".join(x)+"\n")
    f.close()

if __name__ == '__main__':
    makebamwellID(args.sam, args.barcode)
