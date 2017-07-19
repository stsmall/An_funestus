#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 13:04:01 2017
python vcf_fill_gvcf.py -v vcfin
@author: scott
"""
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcf", type=str, required=True,
                    help='merged vcf file')
args = parser.parse_args()


def get_missing(vcfin):
    """Reads in the passed VCF file and records positions with "./.""

    Parameters
    ------
    vcf: vcffile

    Returns
    ------
    missdict: dictionary, dictionary of missing sites
    """

    missdict = defaultdict(lambda: defaultdict(list))
    with open(vcfin, 'r') as vcfin:
        for line in vcfin:
            if not line.startswith("##"):
                if line.startwith("#CHROM"):
                    pop_iix = line.strip().split()[9:]
                else:
                    x = line.strip.split()
                    chrom = x[0]
                    pos = int(x[1])
                    miss = ["./.:" in i for i in x]
                    [missdict[pop_iix[i]][chrom].append(pos) for i,
                     p in enumerate(miss) if p]
    return(missdict, pop_iix)


def get_gvcf(missdict, pop_iix):
    """Uses missdict to collect sites from gvcf file

    Parameters
    ------
    missdict: dictionary

    Returns
    ------
    gvcfdict: dictionary

    """
    gvcfdict = defaultdict(lambda: defaultdict(dict))
    for sample in missdict.keys():
        with open("{}.g.vcf".format(sample)) as gvcf:
            for line in gvcf:
                if not line.startswith("#"):
                    x = line.strip().split()
                    chrom = x[0]
                    spos = int(x[1])
                    info = x[7]
                    if "END" in info:
                        send = info.split("=")
                        pos = range(spos, send+1)
                        for p in pos:
                            if p in missdict[sample][chrom]:
                                gt = x[pop_iix.index(sample) + 9]
                                gvcfdict[sample][chrom].append((pos, gt))
                    else:
                        pos = spos
                        if pos in missdict[sample][chrom]:
                            fill = pop_iix
                            gvcfdict[sample][chrom].append((pos, fill))
    return(gvcfdict)


def make_vcf(vcf, missdict, gvcfdict, pop_iix):
    """Takes as input the missdict since it has positions and the gvcfdict
    since it has the missing fields, and makes a new vcf

    Parameters
    ------
    vcf: file
    missdict: dict, from get_missing
    gvcfdict: dict, from get_gvcf

    Returns
    ------
    vcf: file, filled vcf

    """
    chrom = ''
    with open("{}.filled".format(vcf), 'w') as fill:
        with open(vcf, 'r') as vcfin:
            for line in vcfin:
                if line.startswith("##"):
                    fill.write(line)
                elif line.startswith("#CHROM"):
                    fill.write(line)
                else:
                    x = line.strip().split()
                    if chrom != x[0]:
                        ml = [missdict[key][chrom] for key in missdict.keys()]
                        ml2 = [j for i in ml for j in i]
                        mlist = sorted(set(ml2))
                    chrom = x[0]
                    pos = x[1]
                    if pos in mlist:
                        # find all missing index
                        miss_iix = [i for i, j in enumerate(x) if "./." in j]
                        for m in miss_iix:
                            sample = pop_iix[m - 9]
                            gtfill = gvcfdict[sample][chrom][pos]
                            x[m] = gtfill
                        fill.write("{}\n".format("\t".join(x)))
                    else:
                        fill.write(line)
    return(None)

if __name__ == "__main__":
    vcfin = args.vcf
    md, pop = get_missing(vcfin)
    gd = get_gvcf(md, pop)
    make_vcf(vcfin, md, gd, pop)
