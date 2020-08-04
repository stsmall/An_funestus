#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:20:03 2020
@author: Scott T. Small

This program calculates a single summary line of mean or median from a file of
pairwise statistics.

Example
-------

   $ for c in *.FK.noncoding.stats;do  # pop size 12 and 10
         grep -v "inf" ${c} | grep -v "nan" > ${c}.2
         python filetObs.py -i ${c}.2 --h1 12 --h2 10 --mask >> stats.mask.out2
     done

Notes
-----

    It was designed specifically to work on the out file from FILET analysis.
    The obs function inf abc_stats will return a single vector as well
    so this program is not needed.

"""

import sys
import argparse
import numpy as np


def getStats(inFile, h1, h2, mask):
    """Check and summarize a file of window-based stats for a pairwise comparison.

    Parameters
    ----------
    inFile : TYPE
        DESCRIPTION.
    h1 : TYPE
        DESCRIPTION.
    h2 : TYPE
        DESCRIPTION.
    mask : TYPE
        DESCRIPTION.

    Returns
    -------
    header : TYPE
        DESCRIPTION.
    stats_list : TYPE
        DESCRIPTION.

    """
    stats_list = []
    header = []
    keep_stats = np.array([True, True, True, True, False, True, False, False,
                           True, True, True, True, True, False, True, False,
                           False, True, True, True, True, True, True, True,
                           True, True, True, True, False, False, False])
    with open(inFile, 'r') as stats:
        for line in stats:
            if line.startswith("chrom"):
                header = np.array(next(stats).split()[4:])
            elif line.startswith("pi"):
                header = np.array(next(stats).split())
            else:
                x = line.split()
                stat_arr = np.array(x[4:], dtype=np.float)
                if not any(np.isinf(stat_arr)):
                    filt = (stat_arr <= 1)
                    filt_1 = np.ones(len(stat_arr))
                    # tajD
                    if not filt[5]:
                        if -5 < stat_arr[5] < 5:
                            filt[5] = True
                    if not filt[14]:
                        if -5 < stat_arr[14] < 5:
                            filt[14] = True
                    # hap counts
                    if stat_arr[7] <= h1:
                        filt[7] = True
                    if stat_arr[16] <= h2:
                        filt[16] = True
                    # zx
                    if not filt[23]:
                        if stat_arr[23] < 5:
                            filt[23] = True
                    # dd1 dd2 dxy_min
                    if filt[21]:
                        filt[24] = True
                        filt[25] = True
                    xfilt = filt * filt_1
                    xfilt[xfilt == 0] = np.nan
                    stats = stat_arr * xfilt
                    if mask:
                        stats_list.append(stats[keep_stats])
                    else:
                        stats_list.append(stats)
    if mask and header:
        header = header[keep_stats]
    return header, stats_list


def sumStats(header, stats_list, mean):
    if mean:
        m = np.nanmean(np.vstack(stats_list), axis=0)
    else:
        m = np.nanmedian(np.vstack(stats_list), axis=0)
    print(f"{' '.join(map(str, header))}\n{' '.join(map(str, m))}")


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', "--inFile", required=True, help="tab delimted file "
                        "expects the first 4 columns to be Chrom\tStart\tEnd\tSites")
    parser.add_argument('--h1', type=int, required=True, help="sample size in haps "
                        "of first population")
    parser.add_argument('--h2', type=int, required=True, help="sample size in haps "
                        "of second population")
    parser.add_argument("--mean", action="store_true", help="returns the mean instead "
                        "of the median")
    parser.add_argument("--mask", action="store_true", help="if added then removes "
                        "statstistics that are haplotype dependent. expect 22 "
                        "instead of 31")
    return parser.parse_args(args_in)


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    infile = args.inFile
    mask = args.mask
    mean = args.mean
    h1 = args.h1
    h2 = args.h2
    # =========================================================================
    #  Main executions
    # =========================================================================
    header, stats_list = getStats(infile, h1, h2, mask)
    sumStats(header, stats_list, mean)


if __name__ == "__main__":
    main()
