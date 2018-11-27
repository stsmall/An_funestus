#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 17:14:44 2018
accessibility.py --bed
@author: scott
"""

from __future__ import print_function
from __future__ import division
import argparse
import glob
from collections import defaultdict
from collections import Counter
import numpy as np
import pysam
import gzip
from math import sqrt
import logging

#logging.debug('This message should appear on the console')
#logging.info('So should this')
#logging.warning('And this, too')

parser = argparse.ArgumentParser()
parser.add_argument('-c', "--chroms", nargs='+', action="append",
                    required=True, help="list of chromosomes")
parser.add_argument("--loglevel", help="info, debug, warning, error, critical")
parser.add_argument("--log", action="store_false")
args = parser.parse_args()


def modeCount(modearray, chrom, covdict):
    """
    """
    # use mode as in Ag1000g
    try:
        mode = Counter(modearray).most_common(1)[0][0]
    except IndexError:
        import ipdb;ipdb.set_trace()
    for pos in np.where(modearray == 0)[0]:  # zeroCov
        try:
            covdict[chrom][pos+1][0] += 1
        except KeyError:
            covdict[chrom][pos+1] = [1, 0, 0, 0, 0]
    for pos in np.where(modearray > mode*2)[0]:  # highCov
        try:
            covdict[chrom][pos+1][1] += 1
        except KeyError:
            covdict[chrom][pos+1] = [0, 1, 0, 0, 0]
    for pos in np.where(modearray < mode/2)[0]:  # lowCov
        try:
            covdict[chrom][pos+1][2] += 1
        except KeyError:
            covdict[chrom][pos+1] = [0, 0, 1, 0, 0]
    return(covdict)


def avgCount(modearray, chrom, covdict):
    """
    """
    # use avg from Li 2012??
    avg = np.mean(modearray)
    for pos in np.where(modearray == 0)[0]:  # zeroCov
        try:
            covdict[chrom][pos][0] += 1
        except KeyError:
            covdict[chrom][pos] = [1, 0, 0, 0, 0]
    for pos in np.where(modearray > 3*sqrt(avg) + avg)[0]:  # highCov
        try:
            covdict[chrom][pos][1] += 1
        except KeyError:
            covdict[chrom][pos] = [0, 1, 0, 0, 0]
    for pos in np.where(modearray < 3*sqrt(avg) + avg)[0]:  # lowCov
        try:
            covdict[chrom][pos][2] += 1
        except KeyError:
            covdict[chrom][pos] = [0, 0, 1, 0, 0]
    return(covdict)


def maskCov(bedlist, chromlist, modefx=True):
    """Percentage of Individuals with coverage min and max at the position

    """
    covdict = defaultdict(dict)
    first_chrom = chromlist[0]
    chrom = first_chrom
    chrlendict = {}
    pos = 1
    logging.info('{}\n'.format(bedlist))
    for bed in bedlist:
        indvbed = gzip.open(bed, 'r')
        modelist = []
        for line in indvbed.readlines():
            x = line.split()
            if line == "":
                chrlendict[chrom] = pos
                modearray = np.array(modelist)
                if modefx:
                    covdict = modeCount(modearray, chrom, covdict)
                else:
                    covdict = avgCount(modearray, chrom, covdict)
            else:
                if chrom == x[0]:  # same chrom
                    modelist.append(int(x[2]))
                else:  # different chrom
                    chrlendict[chrom] = pos
                    modearray = np.array(modelist)
                    if modefx:
                        covdict = modeCount(modearray, chrom, covdict)
                    else:
                        covdict = avgCount(modearray, chrom, covdict)
                    modelist = []
                    modelist.append(int(x[2]))
                chrom = x[0]  # new chrom
                pos = int(x[1])
    return(covdict, chrlendict)


def maskQual(bamlist, chromlendict, covdict):
    """
    """
    step = 10000
    for bam in bamlist:
        bamfile = pysam.AlignmentFile(bam, "rb")
        for chrm in chromlendict.keys():
            start = 1
            stop = step
            while start < chromlendict[chrm]:
                try:
                    cols = bamfile.pileup(chrm, start, stop)
                except IndexError:
                    cols = bamfile.pileup(chrm, start, chromlendict[chrm])
                for piles in cols:
                    if np.mean(piles.get_mapping_qualities()) < 30:  # low mapQual
                        covdict[chrm][piles.reference_pos + 1][3] += 1
                    if sum(piles.get_mapping_qualities()) == 0:  # ambig read, MQ of 0
                        covdict[chrm][piles.reference_pos + 1][4] += 1
                start = stop
                stop += step
        bamfile.close()
    return(covdict)


if __name__ == '__main__':
    # let's try logging for a change
    loglevel = args.loglevel
    if loglevel:
        numeric_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: %s' % loglevel)
        logging.basicConfig(filename='accessibility.log', level=numeric_level,
                            format='%(asctime)s %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.Logger.disabled = args.log
    bedlist = glob.glob("*.genCov.gz")
    chromlist = args.chroms
    covdict, chromlendict = maskCov(bedlist, chromlist[0])
#    bamlist = glob.glob("*mdup.bam")
#    covdict = maskQual(bamlist, chromlendict, covdict)
    with open("accessiblity.pos.txt", 'w') as posout:
        for k in covdict.keys():
            for p in covdict[k]:
                posout.write("{}\t{}\t{}\n".format(k, p, "\t".join(map(str, covdict[k][p]))))
