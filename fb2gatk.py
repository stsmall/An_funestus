#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 12:07:50 2017
Convert Freebayes VCF format fields to GATK format fields for merging
@author: stsmall
"""
from __future__ import print_function
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, action="store",
                    required=True, help="name of infile")
args = parser.parse_args()


def vcfformat(geno):
    """
    GT:GQ:DP:AD:RO:QR:AO:QA:GL 0/1:99:40:19,11:19:622:11:399:-2.8296,0,-24.7889

    GT:AD:DP:GQ:PL 0/0:994,0:994:99:0,120,1800
    """
    try:
        gt = geno.split(":")
        # transfer is fb format, ref is GATK
        ad = gt[3]
        dp = gt[2]
        gq = int(float(gt[1]))
        gl = gt[-1].split(",")  # from GL to PL
        if "." in gt[0]:
            gt = "{}:{}:{}:{}:{}".format(gt[0], ad, dp, gq, gt[-1])
        else:
            raw_pl = [-10 * float(i) for i in gl]
            norm_pl = min(raw_pl)
            pl = [int(i - norm_pl) for i in raw_pl]
            plstr = map(str, pl)
            gt = "{}:{}:{}:{}:{}".format(gt[0], ad, dp, gq, ",".join(plstr))
    except ValueError:
        import ipdb;ipdb.set_trace()
    return(gt)


def loadvcf(vcfFile):
    """
    """
    outstream = open("{}.gatk".format(vcfFile), 'w')
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                outstream.write(line)
            else:
                x = line.strip().split()
                for sample in range(9, len(x)):
                    geno = vcfformat(x[sample])
                    x[sample] = geno
                x[8] = "GT:AD:DP:GQ:PL"
                outstream.write("{}\n".format("\t".join(x)))
    outstream.close()
    return(None)


if __name__ == "__main__":
    loadvcf(args.vcfFile)
