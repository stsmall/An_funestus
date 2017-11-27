#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 13:28:53 2017

@author: scott
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('INvcf', metavar="INvcf", type=str,
                    help='path to vcf file')
args = parser.parse_args()


def fixPGTPID(vcf):
    """the PGT and PID fields are not always in the same order
       this function removes them since they are a pain in the ass. if they
       contain phase information, that info is copied to the genotype
    """
    f = open(vcf + '.fix', 'w')
    with open(vcf, 'r') as vcffile:
        for line in vcffile:
            if line.startswith("##") or line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                if "PGT" in x[8]:
                    if "PL:PGT:PID" in x[8]:
                        for sample in range(9, len(x)):
                            gt = x[sample].split(":")
                            newgt = [gt[0], gt[1], gt[2], gt[3], gt[4]]
                            x[sample] = ":".join(newgt)
                        x[8] = "GT:AD:DP:GQ:PL"
                        f.write("{}\n".format("\t".join(x)))
                    elif "PGT:PID:PL" in x[8]:
                        for sample in range(9, len(x)):
                            gt = x[sample].split(":")
                            newgt = [gt[0], gt[1], gt[2], gt[3], gt[6]]
                            x[sample] = ":".join(newgt)
                        x[8] = "GT:AD:DP:GQ:PL"
                        f.write("{}\n".format("\t".join(x)))
                elif "." in x[4]:
                    # fix invariant
                    # GT:AD:DP:RGQ from HaplotypeCaller
                    for sample in range(9, len(x)):
                        gt = x[sample].split(":")
                        newgt = [gt[0], gt[1], gt[2], gt[3], "0,500,500"]
                        x[sample] = ":".join(newgt)
                    f.write("{}\n".format("\t".join(x)))
                else:
                    f.write(line)
    f.close()
    return(None)


if __name__ == '__main__':
    fixPGTPID(args.INvcf)
