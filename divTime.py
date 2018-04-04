#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 17:41:02 2018
Uses method from Theunert and Slatkin : Estimation of population divergence
times from SNP data and a test for treeness.doi: https://doi.org/10.1101/281881

Usage:
divTime.py -v VCF -p 1 2 -p 5 6 -o 10 -msmc FILE [--mle] [--boots INT]

IMPORTANT: the file passed as -msmc should be in an ms formatted line.
You should use MSMC2ms-sts.py --msmc FOO > ms_format.out to make this file.
script is available at https://github.com/stsmall/Wb_sWGA

@author: stsmall
"""
import numpy as np
import allel
from math import log as log
from math import exp as exp
from scipy.optimize import minimize
from itertools import combinations
import os.path
import h5py
import argparse
import sys
if sys.version_info[0] < 3:
    raise "Must be using Python 3"

# TODO: MLE check on k0, k1, c; email Slatkin

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, required=True,
                    help="vcf file of variants")
parser.add_argument('-p', "--pops", nargs='+', action="append",
                    required=True, help="index of populations"
                    "-p 1 2 3 -p 6 7 8")
parser.add_argument('-o', "--outgroup", nargs='+', action="append",
                    help="index of outgroup")
parser.add_argument('-msmc', "--piecewise", help="use psmc for population"
                    "size changes, if None then assumes constant. First run"
                    "MSMC2ms-sts.py --msmc FOO")
parser.add_argument("--boots", type=int, default=0,
                    help="number of bootstraps")
parser.add_argument("--mle", action="store_true", help="use MLE, default is"
                    "fast approx. MLE requires larger popsizes")
args = parser.parse_args()


def loadvcf(vcfFile):
    """Reads VCF using scikit-allel, object is stored as pandas DF
    """
    print("loading vcf file...")
    print("using scitkit allele version:", allel.__version__)
    h5 = "{}h5".format(vcfFile.strip("vcf"))
    if os.path.isfile(h5):
        callset = h5py.File(h5, 'r')
    else:
        callset = allel.read_vcf(vcfFile)
        print("creating h5 for faster loading")
        allel.vcf_to_hdf5(vcfFile, h5)
    return(callset)


def filterGT(callset, pops, outgroup):
    """Count patterns from VCF
    """
    gt = allel.GenotypeArray(callset['calldata/GT'])
    if outgroup:
        # filter on outgroup pop
        acs = gt[:, outgroup].count_alleles(max_allele=1)
        flt = acs.is_segregating()
    else:
        # filter without using outgroup using sampled pops
        subpops = {"popA": pops[0],
                   "popB": pops[1]
                   }
        acs = gt.count_alleles_subpops(subpops, max_allele=1)
        acu = allel.AlleleCountsArray(acs["popA"][:] + acs["popB"][:])
        flt = acu.is_segregating()
    # remove non-segrating
    gt = gt.compress(flt, axis=0)
    return(gt)


def estimCandK(n1, n2, n3, mle):
    """ML estimates of k and c
    """
    k_hat = .75 * ((2*n3 + n2) / (n1 + n2 + n3))
    c_hat = (2*n3 - n2) / (2*n3 + n2)
    if mle:
        fun = lambda x: -n1*np.log((1-x[1])+(x[0]*x[1])/3) + -n2*np.log((2*(1-x[0])*x[1])/3) + -n3*np.log(((1+x[0])*x[1])/3)
        bnds = ((1E-9, 1-1E-9), (1E-9, 1-1E-9))
        cons = [{'type': 'ineq', 'fun': lambda x: (1-x[1])+((x[0]*x[1])/3) - 1E-9},
                {'type': 'ineq', 'fun': lambda x: (2*(1-x[0])*x[1])/3 - 1E-9},
                {'type': 'ineq', 'fun': lambda x: ((1+x[0])*x[1])/3 - 1E-9}]
        res = minimize(fun, (c_hat, k_hat), method="SLSQP", bounds=bnds)
        c, k = res.x
    else:
        c = c_hat
        k = k_hat
    return(c, k)


def countN1N2N3(gt, pops, mle):
    """
    """
    clist = []
    klist = []
    # make gt arrays for each subpop, then haplotype arrays
    gtA = gt.take(pops[0], axis=1)
    htA = gtA.to_haplotypes()
    gtB = gt.take(pops[1], axis=1)
    htB = gtB.to_haplotypes()
    if len(pops[1]) == 1:
        hap2 = list(range(len(pops[1]*2)))
        for hap1 in list(range(len(pops[0]))):
            ma = htA[:, [hap1]].count_alleles(max_allele=1)
            mb = htB[:, hap2].count_alleles(max_allele=1)
            jsfs = allel.joint_sfs(ma[:, 1], mb[:, 1])
            try:
                n1 = jsfs[0, 2] + jsfs[1, 0]
                n2 = jsfs[0, 1] + jsfs[1, 1]
                n3 = jsfs[0, 0] + jsfs[1, 2]
            except IndexError:
                z = np.zeros((2, 3), dtype=int)
                z[:jsfs.shape[0], :jsfs.shape[1]] = jsfs
                n1 = z[0, 2] + z[1, 0]
                n2 = z[0, 1] + z[1, 1]
                n3 = z[0, 0] + z[1, 2]
            c_hat, k_hat = estimCandK(n1, n2, n3, mle)
            clist.append(c_hat)
            klist.append(k_hat)
    else:
        for hap1 in list(range(len(pops[0]))):
            for hap2 in list(combinations(range(len(pops[1])*2), 2)):
                ma = htA[:, [hap1]].count_alleles(max_allele=1)
                mb = htB[:, hap2].count_alleles(max_allele=1)
                z = np.zeros((2, 3), dtype=int)
                jsfs = allel.joint_sfs(ma[:, 1], mb[:, 1])
                try:
                    n1 = jsfs[0, 2] + jsfs[1, 0]
                    n2 = jsfs[0, 1] + jsfs[1, 1]
                    n3 = jsfs[0, 0] + jsfs[1, 2]
                except IndexError:
                    z = np.zeros((2, 3), dtype=int)
                    z[:jsfs.shape[0], :jsfs.shape[1]] = jsfs
                    n1 = z[0, 2] + z[1, 0]
                    n2 = z[0, 1] + z[1, 1]
                    n3 = z[0, 0] + z[1, 2]
                c_hat, k_hat = estimCandK(n1, n2, n3, mle)
                clist.append(c_hat)
                klist.append(k_hat)
    return(np.mean(clist), np.mean(klist))


def estimDiv(c, psmc, r, t):
    """Estimate divergence using eq 12
    """
    N0 = 0
    if psmc:
        if not r:
            # parse psmc
            f = open(psmc, 'r')
            line = f.readline().split("-eN ")
            t = [float(i.split()[0]) for i in line[1:]]
            t.insert(0, 0.0)
            r = [float(i.split()[1]) for i in line[1:]]
            N0 = float(line[0].split()[1]) / float(line[0].split()[4])
            r.insert(0, 1.0)
        i = 0
        nc = 1.0
        while (1-nc*exp(-(t[i+1]-t[i])/r[i])) < c:
            nc *= exp(-(t[i+1]-t[i])/r[i])
            i += 1
            #print("i:{}, t[i]:{}, t[i+1]:{}, r[i]:{}, nc:{}".format(i, t[i], t[i+1], r[i], nc))
        j = i
        print("nc = {}, 1-nc = {}".format(nc, 1-nc))
        T_hat = -r[j]*log((1-c) / nc) + t[j]
    else:
        T_hat = -log(1-c)  # assumes constant popsize
    return(r, t, N0, T_hat)


def calcCI(gt, pops, psmc, boots, r, t):
    """
    """
    print("Running bootstraps...")
    T_hatlist = []
    # random resampling
    indices_rs = np.random.randint(0, len(gt), (1, boots, len(gt)))
    for b in range(boots):
        print("bootstrap number {}".format(b+1))
        gt = gt.take(indices_rs[0][b], axis=0)
        c, k = countN1N2N3(gt, pops)
        r, t, N0, T_hat = estimDiv(c, psmc, r, t)
        T_hatlist.append(T_hat)
    # quantiles
    t_lowCI = np.percentile(T_hatlist, 0.025)
    t_highCI = np.percentile(T_hatlist, 0.975)
    return(t_lowCI, t_highCI)


if __name__ == "__main__":
    popset = args.pops
    msmc = args.piecewise
    # pop_ix = [tuple(map(int, x)) for x in popset]
    pop_ix = [list(map(int, x)) for x in popset]
    if args.outgroup:
        outgroup_ix = [list(map(int, x)) for x in args.outgroup][0]
    else:
        outgroup_ix = ''
    callset = loadvcf(args.vcfFile)
    gt = filterGT(callset, pop_ix, outgroup_ix)
    c, k = countN1N2N3(gt, pop_ix, args.mle)
    print("c_hat = {}; k1_hat = {}".format(c, k))
    r = ''
    t = ''
    r, t, N0, T_hat = estimDiv(c, msmc, r, t)
    if args.boots > 0:
        t_LCI, t_HCI = calcCI(gt, pop_ix, msmc, args.boots, r, t)
        print("{} in 2Ne gens ({} - {})".format(T_hat, t_LCI, t_HCI))
    else:
        print("{} in 2Ne gens".format(T_hat))
    if N0 > 0:
        print("theta at time 0: {}".format(N0))
