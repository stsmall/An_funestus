#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:20:03 2020
@author: Scott T. Small

"""
import sys
import argparse
import numpy as np


def getStats(inFile, h1, h2, mask):
    stats_list = []
    keep_stats = np.array([True, True, True, True, False, True, False, False,
                           True, True, True, True, True, False, True, False,
                           False, True, True, True, True, True, True, True,
                           True, True, True, True, False, False, False])
    with open(inFile, 'r') as stats:
        header = np.array(next(stats).split()[4:])
        for line in stats:
            x = line.split()
            stat_arr = np.array(x[4:], dtype=np.float)
            if not any(np.isinf(stat_arr)):
                breakpoint()
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
                        filt[5] = True
                xfilt = filt * filt_1
                xfilt[xfilt == 0] = np.nan
                stats = stat_arr * xfilt
                if mask:
                    stats_list.append(stats[keep_stats])
                else:
                    stats_list.append(stats)
    if mask:
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
    parser.add_argument('-i', "--inFile", required=True)
    parser.add_argument('--h1', type=int, required=True)
    parser.add_argument('--h2', type=int, required=True)
    parser.add_argument("--mean", action="store_true")
    parser.add_argument("--mask", action="store_true")
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
