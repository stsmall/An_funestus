#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 15:17:54 2017
fillrefinfb.py
@author: stsmall
"""
import argparse
import ipdb
parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf file')
parser.add_argument('-o', "--outgroup", type=str,
                    help='outgroup name')
args = parser.parse_args()


def fixref(vcfFile, outgroup):
    """
    """
    f = open("{}.reffill".format(vcfFile), 'w')
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if line.startswith("##"):
                f.write(line)
            elif line.startswith("#CHROM"):
                samplelist = line.strip().split()
                x = line.strip().split()
                sample_ix = range(9, len(x))
                outgroups_ix = [i for i, y in enumerate(samplelist) if outgroup in y.split(".")[0]]
                ingroups_ix = [var for var in sample_ix if var not in outgroups_ix][9:]
                f.write(line)
            else:
                x = line.strip().split()
                gtout = [x[i].split(":")[0] for i in outgroups_ix]
                if gtout.count("./.") < len(outgroups_ix):
                    for s in ingroups_ix:
                        gt = x[s].split(":")
                        if "./." or "." in gt[0]:
                            gt[0] = "0/0"
                            x[s] = ":".join(gt)
                    f.write("{}\n".format("\t".join(x)))
                else:
                    f.write(line)
    f.close()
    return(None)


if __name__ == '__main__':
    fixref(args.INvcf, args.outgroup)
