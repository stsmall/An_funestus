#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 11:35:48 2017
fasta2vcf.py
@author: stsmall
"""
from __future__ import print_function
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', "--fastaFile", type=str,
                    required=True, help="name of infile")
args = parser.parse_args()

vcfdict = {}
inds = ['P1', 'P2', 'P3', 'P4', 'PO']
with open(args.fastaFile, 'r') as fasta:
    for line in fasta:
        if line.startswith(">"):
            header = line.strip().lstrip(">")
            seq = fasta.next()
            vcfdict[header] = seq
            sites = len(seq)
f = open('vcf.out', 'w')
pos = 0
for seq in range(sites):
    outgroup = vcfdict['PO'][seq]
    gt = []
    for sample in inds:
        base = vcfdict[sample][seq]
        if base == outgroup:
            gt.append("0/0:.:.:.")
        else:
            gt.append("1/1:.:.:.")
    geno = "\t".join(gt)
    f.write("X\t{}\t.\tA\tT\t50\tPASS\t.\t.\t{}\n".format(pos, geno))
    pos += 1
f.close()
