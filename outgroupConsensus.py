#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 16:21:43 2018

@author: scott
"""

import sys
import argparse
assert sys.version_info < (3, 0)

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf',
                    help='path to vcf', required=True)
parser.add_argument('-o', "--outgroup", required=True,
                    help="index of outgroup")
args = parser.parse_args()


def collapseOutgroup(vcfFile, outgroup_ix):
    """Combine the gt calls of outgroup individual is fixed and 1 is missing
    """
    f = open("{}.outgroup".format(vcfFile), 'w')
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if line.startswith("##"):
                f.write(line)
            elif line.startswith("#CHROM"):
                samples = line.split()
                samples.append("outgroup")
                f.write("{}\n".format('\t'.join(samples)))
            else:
                x = line.split()
                gt = []
                for ix in outgroup_ix:
                    gt.append(x[ix+9])
                # do something with gt
                gt_out = "\t".join(gt)
                homR = gt_out.count("0/0")
                homA = gt_out.count("1/1")
                het = gt_out.count("0/1")
                if homR + homA + het == 0:
                    gt_con = "./.:.:.:.:."
                else:
                    o_ix = [homR, homA, het].index(max([homR, homA, het]))
                    if o_ix == 0:
                        gt_con = next(x for x in gt if '0/0' in x)
                    elif o_ix == 1:
                        gt_con = next(x for x in gt if '1/1' in x)
                    else:
                        gt_con = next(x for x in gt if '0/1' in x)
                x.append(gt_con)
                f.write("{}\n".format('\t'.join(x)))
    return(None)


if __name__ == "__main__":
    out_ix = map(int, args.outgroup)
    collapseOutgroup(args.vcf, out_ix)
