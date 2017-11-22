#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 14:04:16 2017
@author: stsmall
geno2vcf
"""
import argparse
from collections import defaultdict
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--vcf", type=str, required=True,
                    help="infile in vcf")
parser.add_argument('-g', "--geno", type=str, required=True,
                    help="infile in geno")
parser.add_argument('-f', "--fasta", type=str, required=True,
                    help="bed file from a fasta file, 1 base per site. Make"
                    "this by first making a bed file of all positions and then"
                    "using bedtools getFastafrombed -tab -Outbed")
args = parser.parse_args()


def readbedfasta(fasta):
    """
    """
    fdict = defaultdict(dict)
    with open(fasta, 'r') as fbed:
        for line in fbed:
            x = line.strip().split()
            fdict[x[0]][x[2]] = x[-1]
    return(fdict)


def getvcfheaderinfo(vcf):
    """
    """
    vcfdict = defaultdict(dict)
    header = []
    with open(vcf, 'r') as v:
        for line in v:
            if line.startswith("##"):
                header.append(line)
            elif line.startswith("#CHROM"):
                samples = line.strip().split()
            else:
                x = line.strip().split()
                vcfdict[x[0]][x[1]] = [x[3], x[4], x[5]]
    return(header, samples, vcfdict)


def geno2vcf(geno, vcfdict, fdict, header, samples):
    """
    """
    with open("{}.vcf".format(geno), 'w') as f:
        for h in header:
            f.write(h)
        f.write("{}\n".format("\t".join(samples)))
        with open(geno, 'r') as geno:
            for line in geno:
                if not line.startswith("#"):
                    x = line.strip().split()
                    chrom = x[0]
                    pos = x[1]
                    try:
                        ref = vcfdict[chrom][pos][0]
                        alt = vcfdict[chrom][pos][1]
                        qual = vcfdict[chrom][pos][2]
                    except KeyError:
                        ref = fdict[chrom][pos]
                        alt = ''
                        for r in x[1:]:
                            if r != "{}/{}".format(ref, ref):
                                if "N" not in r:
                                    alt = r.split("/")[0]
                        if alt == '':
                            alt = 'N'
                        qual = '100'

        f.write("{}\t{}\t.\t{}\t{}\t{}\tPASS\t.\tGT\t{}\n".format(chrom, pos, ref, alt, qual, "\t".join(genos)))









