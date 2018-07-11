#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 13:44:12 2018

@author: scott
"""

from __future__ import print_function
from __future__ import division
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, required=True,
                    help="vcf file")
parser.add_argument('-e', "--estFile", type=str, required=True,
                    help="est-sfs output with 1st and 2nd columns having CHR"
                    "POS")
args = parser.parse_args()


def readEstSFS(estFile):
    """Loads estFile in dictionary
    """
    estdict = {}
    with open(estFile, 'r') as est:
        for line in est:
            x = line.split()
            chrpos = "{}_{}".format(x[0], x[1])
            # probMajor probA probC probG probT
            estdict[chrpos] = x[4:]
    return(estdict)


def derivedVCF(estdict, vcfFile):
    """
    """
    nucstr = 'ACGT'
    a = open("ancAllel.out", 'w')
    f = open("{}.derived".format(vcfFile), 'w')
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                ref = x[3]
                alt = x[4]
                # gt only and count major
                calt = 0
                cref = 0
                for i, sample in enumerate(x[8:]):
                    gt = sample.split(":")[0]
                    calt += gt.count('1')
                    cref += gt.count('0')
                    x[i+8] = gt
                if cref > calt:
                    maj = ref
                else:
                    maj = alt
                # find derived
                chrpos = "{}_{}".format(x[0], x[1])
                estlist = map(float, estdict[chrpos])
                # estlist = [pMaj, pA, pC, pG, pT]
                if estlist[0] >= 0.70:
                    if maj == ref:
                        f.write("{}\n".format("\t".join(x)))
                    elif maj == alt:
                        x[4] = ref
                        x[3] = alt
                        for i, gt in enumerate(x[8:]):
                            if gt == '0/0':
                                x[i+8] = '1/1'
                            elif gt == '1/1':
                                x[i+8] = '0/0'
                        f.write("{}\n".format("\t".join(x)))
                        a.write("{}\t{}\t{}\n".format(x[0], x[1], maj))
                else:
                    pNuc = max(estlist[1:])
                    if pNuc >= 0.70:
                        ix = estlist[1:].index(pNuc)
                        nuc = nucstr[ix]
                        if nuc == ref:
                            f.write("{}\n".format("\t".join(x)))
                        elif nuc == alt:
                            x[4] = ref
                            x[3] = alt
                            for i, gt in enumerate(x[8:]):
                                if gt == '0/0':
                                    x[i+8] = '1/1'
                                elif gt == '1/1':
                                    x[i+8] = '0/0'
                            f.write("{}\n".format("\t".join(x)))
                            a.write("{}\t{}\t{}\n".format(x[0], x[1], nuc))
                    else:
                        pass
    f.close()
    a.close()
    return(None)


if __name__ == "__main__":
    estdict = readEstSFS(args.estFile)
    derivedVCF(estdict, args.vcfFile)
