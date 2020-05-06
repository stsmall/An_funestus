#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 15:08:11 2018

##awk 'OFS = "\t",$4 = $3 + 1' BEDGRAPH | cut -f 1,3- > QUERY.bed
sed '1d' BEDGRAPH > BED
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
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(x[10], x[9][0],
                                                              x[12], x[13],
                                                              x[14], x[9][1],
                                                              x[16], x[17]))
    f.close()
    return(None)


if __name__ == "__main__":
    psl2tln(args.INpsl)
