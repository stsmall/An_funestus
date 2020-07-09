#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 13:04:01 2017
python vcf_fill_gvcf.py -v vcfin
@author: scott
"""
from collections import defaultdict
import glob
import argparse
import progressbar as pb
import numpy as np
_widgets = [pb.Percentage(), pb.Bar()]
parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcf", type=str, required=True,
                    help='merged vcf file')
args = parser.parse_args()


def bufcount(filename):
    f = open(filename)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read  # loop optimization
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
    return lines


def get_filtered():
    """
    """
    vcf_files = glob.glob("*.snps.flt.vcf")
    n = len(vcf_files)
    fltdict = defaultdict(lambda: defaultdict(list))
    # progress bar
    progress = pb.ProgressBar(widgets=_widgets, maxval=n).start()
    progvar = 0
    print("Getting Filtered Sites\n")
    for v in vcf_files:
        with open(v, 'r') as vcf:
            for line in vcf:
                if not line.startswith("##"):
                    if line.startswith("#CHROM"):
                        sample = line.strip().split()[9]
                    else:
                        x = line.strip().split()
                        if x[6] != "PASS":
                            chrom = x[0]
                            pos = int(x[1])
                            fltdict[sample][chrom].append(pos)
        progress.update(progvar + 1)
        progvar += 1
    progress.update(n)
    return(fltdict)


def get_missing(vcfin, flt, nlines):
    """Reads in the passed VCF file and records positions with "./.""

    Parameters
    ------
    vcf: vcffile

    Returns
    ------
    missdict: dictionary, dictionary of missing sites
    """
    missdict = defaultdict(lambda: defaultdict(list))
    # progress bar
    progress = pb.ProgressBar(widgets=_widgets, maxval=nlines).start()
    progvar = 0
    print("Getting Missing Sites\n")
    with open(vcfin, 'r') as vcf:
        l = 0
        for line in vcf:
            if not line.startswith("##"):
                if line.startswith("#CHROM"):
                    pop_iix = line.strip().split()[9:]
                else:
                    l += 1
                    x = line.strip().split()
                    chrom = x[0]
                    pos = int(x[1])
                    miss = ["./.:" in i for i in x[9:]]
                    [missdict[pop_iix[i]][chrom].append(pos) for i,
                     p in enumerate(miss) if p]
                    if l % 10000 == 0:
                        # print("done reading line {}\n".format(l))
                        try:
                            progress.update(progvar + 10000)
                            progvar += 10000
                        except ValueError:
                            progress.update(nlines)
    # remove filtered sites
    print("Removing Filtered Sites\n")
    n = len(missdict.keys())
    progress = pb.ProgressBar(widgets=_widgets, maxval=n).start()
    progvar = 0
    for sample in missdict.keys():
        for chrom in missdict[sample].keys():
            a = missdict[sample][chrom]
            b = flt[sample][chrom]
            filtered = list(set(a).symmetric_difference(b))
            missdict[sample][chrom] = np.array(filtered)
        progress.update(progvar + 1)
        progvar += 1
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
    n = len(missdict.keys())
    progress = pb.ProgressBar(widgets=_widgets, maxval=n).start()
    progvar = 0
    print("Getting data from gVCF ... this will take a while\n")
    for sample in missdict.keys():
        with open("{}.g.vcf".format(sample), 'r') as gvcf:
            print("Opening {}.g.vcf".format(sample))
            for line in gvcf:
                if not line.startswith("#"):
                    x = line.strip().split()
                    chrom = x[0]
                    spos = int(x[1])
                    info = x[7]
                    if "END" in info:
                        send = int(info.split("=")[1])
                    else:
                        send = spos
                    a = missdict[sample][chrom]
                    gt_pos = a[np.where(np.logical_and(a >= spos, a <= send))]
                    if gt_pos.any():
                        gt = x[9].split(":")  # GT:DP:GQ:MIN_DP:PL
                        for p in gt_pos:
                            gt2 = "{}:0,{}:{}:{}:{}".format(gt[0], gt[3],
                                                            gt[1], gt[2],
                                                            gt[4])
                            gvcfdict[sample][chrom][str(p)] = gt2
        progress.update(progvar + 1)
        progvar += 1
    return(gvcfdict)


def make_vcf(vcfin, missdict, gvcfdict, pop_iix, nlines):
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
    l = 0
    progress = pb.ProgressBar(widgets=_widgets, maxval=nlines).start()
    progvar = 0
    with open("{}.filled".format(vcfin), 'w') as fill:
        with open(vcfin, 'r') as vcf:
            for line in vcf:
                if line.startswith("##"):
                    fill.write(line)
                elif line.startswith("#CHROM"):
                    fill.write(line)
                else:
                    # progress bar
                    l += 1
                    if l % 10000 == 0:
                        try:
                            progress.update(progvar + 10000)
                            progvar += 10000
                        except ValueError:
                            progress.update(nlines)
                    # fill code
                    x = line.strip().split()
                    if chrom != x[0]:
                        ml = [missdict[key][chrom] for key in missdict.keys()]
                        ml2 = [j for i in ml for j in i]
                        mlist = np.array(sorted(set(ml2)))
                    chrom = x[0]
                    pos = int(x[1])
                    # quick search for pos in array
                    if np.where(mlist == pos)[0]:
                        # find all missing index
                        miss_iix = [i for i, j in enumerate(x[9:])
                                    if "./." in j]
                        for m in miss_iix:
                            sample = pop_iix[m]
                            gtfill = gvcfdict[sample][chrom][str(pos)]
                            x[m + 9] = gtfill
                        fill.write("{}\n".format("\t".join(x)))
                    else:
                        fill.write(line)
    return(None)

if __name__ == "__main__":
    vcfin = args.vcf
    nlines = bufcount(vcfin)
    flt = get_filtered()
    md, pop = get_missing(vcfin, flt, nlines)
    gd = get_gvcf(md, pop)
    make_vcf(vcfin, md, gd, pop, nlines)
