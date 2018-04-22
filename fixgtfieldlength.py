#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 22 15:25:26 2018
Merging from fastphase seems to muck up the field count, causing errors
    during merging. This script simply fixes the gt field length. This func
    should be added to addphased2vcf.py as a check
@author: scott
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf',
                    help='path to vcf', required=True)
args = parser.parse_args()


def fixGTFieldLength(vcfFile):
    """Count the length of the gt field, clip off repeated fields
    """
    f = open("{}.fix".format(vcfFile), 'w')
    with open(vcfFile, 'r') as v:
        for line in v:
            if line.startswith("#"):
                f.write(line)
            else:
                x = line.split()
                for sample in x[9:]:
                    gt = x[sample].split(":")
                    if len(gt) > 5:
                        gt = gt[:5]
                        x[sample] = ":".join(gt)
                    else:
                        pass
                f.write("{}\n".format('\t'.join(x)))
    f.close()
    return(None)


if __name__ == "__main__":
    fixGTFieldLength(args.vcf)
