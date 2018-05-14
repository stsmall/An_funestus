#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 16:21:43 2018

@author: scott
"""

import sys
import argparse
import re
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
    t = open("ancestral_prob.txt", 'w')
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
                gt_out = " ".join(gt)
                homR = len(re.findall(r'0\|0|0/0', gt_out))
                homA = len(re.findall(r'1\|1|1/1', gt_out))
                het = len(re.findall(r'0\|1|0/1', gt_out))
                if (homR + homA + het) == 0:
                    if "|" in line:
                        gt_con = ".|.:.:.:.:."
                    else:
                        gt_con = "./.:.:.:.:."
                else:
                    o_ix = [homR, homA, het].index(max([homR, homA, het]))
                    if o_ix == 0:
                        gt_con = re.search(r'(0\|0|0/0)[^ ]*', gt_out).group()
                        l_ix = x[3]
                    elif o_ix == 1:
                        gt_con = re.search(r'(1\|1|1/1)[^ ]*', gt_out).group()
                        l_ix = x[4]
                    else:
                        gt_con = re.search(r'(0\|1|0/1)[^ ]*', gt_out).group()
                        l_ix = [x[3], x[4]]
                x.append(gt_con)
                lstate = ['A', 'C', 'G', 'T']
                lprob = [0.3, 0.3, 0.3, 0.3]
                if len(l_ix) > 1:
                    for s in l_ix:
                        lprob[lstate.index(l_ix)] = 0.47
                else:
                    lprob[lstate.index(l_ix)] = .93
                t.write("{} {} {}\n".format(x[0], int(x[1])-1, " ".join(map(str, lprob))))
                f.write("{}\n".format('\t'.join(x)))
    t.close()
    f.close()
    return(None)


if __name__ == "__main__":
    out_ix = map(int, args.outgroup)
    collapseOutgroup(args.vcf, out_ix)
