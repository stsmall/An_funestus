#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 14:24:44 2018

python -ms FOO.ms -t 5 --filet /home/pathfilet/

To normalize: cat outfile | python normalizeTwoPopnStats.py None Block > out
	
@author: stsmall
"""
from __future__ import print_function
from __future__ import division
import numpy as np
from subprocess import run, PIPE
import os
from itertools import combinations
import argparse
import sys
import multiprocessing
# check versions
assert sys.version_info[0] >= 3

parser = argparse.ArgumentParser()
parser.add_argument('-ms', "--msFile", type=str, required=True,
                    help="ms formatted file")
parser.add_argument('-t', "--nprocs", type=int, default=os.cpu_count(),
                    help="number of processors")
parser.add_argument("--filet", type=str, help="path to filet exe")
args = parser.parse_args()


def read_msformat_file(msFile):
    """Read and parse simulations from ms formatted file from folder of files
    with basename

    Params
    ------
    msFile: File, msformatted file

    Returns
    ------
    gtlist: list, list of numpy arrays holding haps
    pops: list, list of pop indices
    pos_list: list, list of arrays of positions
    block: int, length of locus

    """
    pos_list = []
    gt_list = []
    with open(msFile, 'r') as ms:
        for line in ms:
            if line != b'':
                # import ipdb;ipdb.set_trace()
                # line = line.decode('utf-8')
                if line.startswith("scrm"):
                    x = line.split()
                    nind = int(x[1])
                    block = int(x[7])
                    pops = int(x[9])
                    popsize = map(int, x[10:10+pops])
                    ind = 0
                    popconfig = []
                    for p in popsize:
                        if p == 0:
                            pass
                        else:
                            popconfig.append(list(range(ind, p+ind)))
                            ind += p
                elif line.startswith("positions"):
                    pos = np.array(line.strip().split()[1:], dtype=np.float64)
                    pos_list.append(pos)
                    # line = next(ms).decode('utf-8')
                    line = next(ms)
                    gt = np.zeros((nind, pos.shape[0]), dtype=np.uint8)
                    cix = 0
                    try:
                        while line:
                            line = list(line.strip())
                            try:
                                gt[cix, :] = np.array(line, dtype=np.uint8)
                            except IndexError:
                                break
                            cix += 1
                            # line = next(ms).decode('utf-8')
                            line = next(ms)
                    except StopIteration:
                        gt_list.append(gt)
                        break
                    gt_list.append(gt)
            else:
                break
    return(gt_list, popconfig, pos_list, block)


def filetStats(gtlist, pops, pos, block, filetpath):
    """Calculate stats from msformatted file

    Params
    ------
    gtlist: list, list of numpy arrays holding haps
    pops: list, list of pop indices
    pos_list: list, list of arrays of positions
    block: int, length of locus
    filetpath: str, path to exe

    Returns
    ------
    None

    """
    loci = len(gtlist)
    for pop1, pop2 in combinations(pops, 2):
        f = open("{}-{}.filetstats.out".format(pop1[0], pop2[0]), 'w')
        n1 = len(pop1)
        n2 = len(pop2)
        fakems = []
        fakems.append("ms {} {} -t tbs -r tbs {} -I 2 {} {}\n1234\n".format(n1+n2, loci, block, n1, n2))
        for i, g in enumerate(gtlist):
            gt = g[pop1+pop2]
            seg_pos = np.sum(gt, axis=0)
            seg_mask = (seg_pos > 0) & (seg_pos < (n1+n2))
            seg = np.count_nonzero(seg_mask)
            posit = pos[i][seg_mask]
            gt_seg = gt[:, seg_mask]
            fakems.append("\n//\nsegsites: {}\npositions: {}\n".format(seg,
                          " ".join(map(str, posit))))
            for a in gt_seg:
                fakems.append("{}\n".format("".join(map(str, a))))
        msinput = "".join(fakems)
        cmd = ["{}twoPopnStats_forML".format(filetpath), str(n1), str(n2)]
        run(cmd, stdout=f, input=msinput, encoding='ascii')
        f.close()
    return(None)


def filetStatsMP(args):
    """Calculate stats from msformatted file
    
    Params
    ------
    args[0]: list, list of pop indices
    args[1]: list, list of pop indices
    args[2]: list, list of numpy arrays holding haps
    args[3]: list, list of arrays of positions
    args[4]: int, length of locus
    args[5]: str, path to exe

    Returns
    ------
    None

    """
    pop1, pop2, gtlist, pos, block, filetpath = args
    loci = len(gtlist)
    with open("{}-{}.filetstats.out".format(pop1[0], pop2[0]), 'w') as f:
        n1 = len(pop1)
        n2 = len(pop2)
        fakems = []
        fakems.append("ms {} {} -t tbs -r tbs {} -I 2 {} {}\n1234\n".format(n1+n2, loci, block, n1, n2))
        for i, g in enumerate(gtlist):
            gt = g[pop1+pop2]
            seg_pos = np.sum(gt, axis=0)
            seg_mask = (seg_pos > 0) & (seg_pos < (n1+n2))
            seg = np.count_nonzero(seg_mask)
            posit = pos[i][seg_mask]
            gt_seg = gt[:, seg_mask]
            fakems.append("\n//\nsegsites: {}\npositions: {}\n".format(seg,
                          " ".join(map(str, posit))))
            for a in gt_seg:
                fakems.append("{}\n".format("".join(map(str, a))))
        msinput = "".join(fakems)
        cmd = ["{}twoPopnStats_forML".format(filetpath), str(n1), str(n2)]
        run(cmd, stdout=f, input=msinput, encoding='ascii')
    return(None)


if __name__ == "__main__":
    gtlist, pops, pos, block = read_msformat_file(args.msFile)
    if args.nprocs == 1:
        filetStats(gtlist, pops, pos, block, args.filet)
    else:
        pool = multiprocessing.Pool(args.nprocs)
        argslist = []
        for pop1, pop2 in combinations(pops, 2):
            argslist.append([pop1, pop2, gtlist, pos, block, args.filet])
        try:
            pool.map(filetStatsMP, argslist)
        except KeyboardInterrupt:
            sys.stdout.write('\033[0m')
            sys.stdout.write('User Interupt\n')
        pool.close()
