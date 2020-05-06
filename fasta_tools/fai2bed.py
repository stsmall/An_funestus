#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 15:37:55 2017
first: python fai2bed.py FASTA.fai
bedtools getfasta -fi FASTA -bed pos.bed -tab -bedOut > FOO.bed
@author: stsmall
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--faiFile", type=str, action="store",
                    required=True, help="name of infile")
args = parser.parse_args()


def makebed(faifile):
    """Make bed coordinates for constructing a bedgraph

    Parameters
    ------
    faifile: file, fasta file index by samtools faidx FASTA

    Returns
    ------
    pos.bed: file, bed file with start/end coordinates

    """
    fdict = {}
    with open(faifile, 'r') as fun:
        for line in fun:
            x = line.strip().split()
            fdict[x[0]] = x[1]
    for chrom in fdict.keys():
        xstart = 0
        f = open("{}.bed".format(chrom), 'w')
        while xstart <= int(fdict[chrom]):
            f.write("{}\t{}\t{}\n".format(chrom, xstart, xstart + 1))
            xstart += 1
        f.close()
    return(None)


if __name__ == "__main__":
    makebed(args.faiFile)
