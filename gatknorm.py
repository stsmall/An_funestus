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
    v = open(vcf + '.trash', 'w')
    f = open(vcf + '.fix', 'w')
    with open(vcf, 'r') as vcffile:
        for line in vcffile:
            if line.startswith("##") or line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                formats = x[8].split(":")
                if x[5] == 'inf':
                    # quality field is INF
                    x[5] = '500'
                if (len(x[3]) > 1):
                    # skip ref allele that are complex insertions
                    v.write(line)
                    pass
                elif any([len(i) > 1 for i in x[4].split(",")]):
                    # skip alt alleles that are complex insertions
                    v.write(line)
                    pass
                elif ("*" not in x[4]) and (len(formats) > 1) and ("<NON_REF>" not in x[4]):
                    # * doesnt actually give the ALT base
                    # *,<NON_REF> does not give ALT base
                    # format fields are ocassionally just GT with no other information
                    if "." in x[4]:
                        # fix invariant
                        for sample in range(9, len(x)):
                            gt = x[sample].split(":")
                            if gt[0] == "./." or gt[0] == ".":
                                x[sample] = "./.:.:.:.:."
                            else:
                                try:
                                    gq = gt[formats.index('RGQ')]
                                except ValueError:
                                    gq = '99'
                                dp = gt[formats.index('DP')]
                                try:
                                    ad = gt[formats.index('AD')]
                                except ValueError:
                                    ad = dp
                                pl = '0'
                                adv = ad.split(",")[0]
                                newgt = [gt[0], adv, dp, gq, pl]
                                x[sample] = ":".join(newgt)
                        x[8] = "GT:AD:DP:GQ:PL"
                        newsite = "\t".join(x)
                        if newsite.count("./.") == len(range(9, len(x))):
                            # every gt is missing
                            v.write(line)
                            pass
                        else:
                            f.write("{}\n".format(newsite))
                    elif "PGT" in formats or "PID" in formats:
                        for sample in range(9, len(x)):
                            gt = x[sample].split(":")
                            if gt[0] == "./." or gt[0] == ".":
                                x[sample] = "./.:.:.:.:."
                            else:
                                ad = gt[formats.index('AD')]
                                dp = gt[formats.index('DP')]
                                gq = gt[formats.index('GQ')]
                                pl = gt[formats.index('PL')]
                                newgt = [gt[0], ad, dp, gq, pl]
                                x[sample] = ":".join(newgt)
                        x[8] = "GT:AD:DP:GQ:PL"
                        newsite = "\t".join(x)
                        if newsite.count("./.") == len(range(9, len(x))):
                            # every gt is missing
                            v.write(line)
                            pass
                        else:
                            f.write("{}\n".format(newsite))
                    else:
                        f.write(line)
                else:
                    v.write(line)
                    pass
    f.close()
    v.close()
    return(None)


if __name__ == '__main__':
    fixPGTPID(args.INvcf)
