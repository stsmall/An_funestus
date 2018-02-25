#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 15:08:11 2018

awk 'OFS = "\t",$4 = $3 + 1' BEDGRAPH | cut -f 1,3- > QUERY.bed
halLiftover --inMemory --outPSLWithName --noDupes FOO.hal QUERY QUERY.bed SUBJECT OUT.bed
PSL2tln.py OUT.bed

@author: stsmall
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('INpsl', metavar="INpsl", type=str,
                    help='path to psl IN file')
args = parser.parse_args()


def psl2tln(pslFile):
    """
    """
    f = open("out.tln", 'w')
    with open(pslFile, 'r') as psl:
        for line in psl:
            x = line.split()
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(x[9], x[8][0],
                                                              x[11], x[12],
                                                              x[13], x[8][1],
                                                              x[15], x[16]))
    f.close()
    return(None)


if __name__ == "__main__":
    psl2tln(args.INpsl)
