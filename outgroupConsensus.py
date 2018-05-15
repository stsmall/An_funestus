#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 16:21:43 2018

@author: scott
"""

import sys
import argparse
import re
import numpy as np
assert sys.version_info < (3, 0)

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf',
                    help='path to vcf', required=True)
parser.add_argument('-o', "--outgroup", nargs='+', action="append",
                    help="index of outgroup")
args = parser.parse_args()


def mutationMatrix(mutarray, anc, ref, alt):
    """Calculates the direction of the mutaion at site
    """
    # this is terribly clumsy, should use something with i !=j for i, j in
    # l = [A, C, G, T]
    if anc == ref:
        if anc == "A":
            if alt == "C":
                mutarray[0, 1] += 1
            elif alt == "G":
                mutarray[0, 2] += 1
            elif alt == "T":
                mutarray[0, 3] += 1
        elif anc == "C":
            if alt == "A":
                mutarray[1, 0] += 1
            elif alt == "G":
                mutarray[1, 2] += 1
            elif alt == "T":
                mutarray[1, 3] += 1
        elif anc == "G":
            if alt == "A":
                mutarray[2, 0] += 1
            elif alt == "C":
                mutarray[2, 1] += 1
            elif alt == "T":
                mutarray[2, 3] += 1
        elif anc == "T":
            if alt == "A":
                mutarray[3, 0] += 1
            elif alt == "C":
                mutarray[3, 1] += 1
            elif alt == "G":
                mutarray[3, 2] += 1
    elif anc == alt:
        if anc == "A":
            if ref == "C":
                mutarray[0, 1] += 1
            elif ref == "G":
                mutarray[0, 2] += 1
            elif ref == "T":
                mutarray[0, 3] += 1
        elif anc == "C":
            if ref == "A":
                mutarray[1, 0] += 1
            elif ref == "G":
                mutarray[1, 2] += 1
            elif ref == "T":
                mutarray[1, 3] += 1
        elif anc == "G":
            if ref == "A":
                mutarray[2, 0] += 1
            elif ref == "C":
                mutarray[2, 1] += 1
            elif ref == "T":
                mutarray[2, 3] += 1
        elif anc == "T":
            if ref == "A":
                mutarray[3, 0] += 1
            elif ref == "C":
                mutarray[3, 1] += 1
            elif ref == "G":
                mutarray[3, 2] += 1
    else:
        pass
    return(mutarray)


def collapseOutgroup(vcfFile, outgroup_ix):
    """Combine the gt calls of outgroup individual is fixed and 1 is missing
    """
    f = open("{}.outgroup".format(vcfFile), 'w')
    t = open("ancestral_prob.txt", 'w')
    mutarray = np.zeros((4, 4))
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
                # import ipdb;ipdb.set_trace()
                for ix in outgroup_ix:
                    gt.append(x[ix+9])
                gt_out = " ".join(gt)
                homR = len(re.findall(r'0\|0|0/0', gt_out))
                homA = len(re.findall(r'1\|1|1/1', gt_out))
                het = len(re.findall(r'0\|1|0/1', gt_out))
                if (homR + homA + het) == 0:
                    l_ix = [x[3], x[4]]
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
                try:
                    lstate = ['A', 'C', 'G', 'T']
                    lprob = [0.03, 0.03, 0.03, 0.03]
                    if len(l_ix) > 1:
                        for s in l_ix:
                            lprob[lstate.index(s)] = 0.47
                    else:
                        lprob[lstate.index(l_ix)] = .93
                        mutarray = mutationMatrix(mutarray, l_ix, x[3], x[4])
                        # well defined AncAllele
                    t.write("{} {} {}\n".format(x[0], int(x[1])-1, " ".join(map(str, lprob))))
                except ValueError:
                    import ipdb;ipdb.set_trace()
                f.write("{}\n".format('\t'.join(x)))
    t.close()
    f.close()
    return(mutarray)


if __name__ == "__main__":
    out_ix = [list(map(int, x)) for x in args.outgroup]
    mutarray = collapseOutgroup(args.vcf, out_ix[0])
    print(mutarray)
