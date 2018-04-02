#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 17:41:02 2018
divTime.py -v VCF -p 1 2 3 4 -p 5 6 7 8 -o 10 11 12 -psmc FILE -N0 10000 --divtime 0.5 [--mle]
uses method from Theunert and Slatkin : Estimation of population divergence
times from SNP data and a test for treeness.doi: https://doi.org/10.1101/281881
@author: stsmall
"""
import numpy as np
import allel
from math import log as log
from math import exp as exp
from scipy.optimize import minimize
from scipy.stats import multivariate_normal as mnorm
from itertools import combinations
import os.path
import h5py
import argparse

# TODO: Numerical solver; PSMC/MSMC conversion
# TODO: MLE check on k0, k1, c; email Slatkin

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, required=True,
                    help="vcf file of variants")
parser.add_argument('-p', "--pops", nargs='+', action="append",
                    required=True, help="index of populations"
                    "-p 1 2 3 -g 6 7 8 -p 11 12 13")
parser.add_argument('-o', "--outgroup", help="index of outgroup")
parser.add_argument('-psmc', "--piecewise", help="use psmc for population"
                    "size changes, if None then assumes constant")
parser.add_argument('-j', "--divtime", help="a boundary where j < Tdiv < j+1"
                    "in coalescent units")
parser.add_argument("--boots", type=int, default=100,
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


def estimation(obs, fun, init, method='Nelder-Mead'):
    mle = lambda param: -np.sum(fun(*[obs, param]))
    result = minimize(mle, init, method=method)
    return result.x


def estimCandK(n1, n2, n3):
    """ML estimates of k and c
    """
    obs = np.array([n1, n2, n3])  # this needs to looklike mvn data
    n3m = np.mean(n3)
    n2m = np.mean(n2)
    n1m = np.mean(n1)
    k_hat = .75 * ((2*n3m + n2m) / (n1m + n2m + n3m))
    c_hat = (2*n3m - n2m) / (2*n3m + n2m)
    init = [k_hat, c_hat]
    # minimize the neg log lk of a mnorm by changing c and k
    # p1 = (1 - k) + ((c*k) / 3)
    # p2 = ((2*(1 - c)*k) / 3)
    # p3 = (((1 + c)*k) / 3)
    c, k = estimation(obs, lambda ob, p: mnorm.logpdf(ob, [(1 - p[0])+((p[1]*p[0]) / 3), ((2*(1 - p[1])*p[0]) / 3), (((1 + p[1])*p[0]) / 3)]), init)
    return(c, k)


def countN1N2N3MLE(gt, pops):
    """Count patterns from VCF
    """
    # use allel to load a vcf, then count patterns of n1,n2,n3 for each
    # combination and estimating k1 and c.
    n1list = []
    n2list = []
    n3list = []
    gtA = gt.take(pops[0], axis=1)
    htA = gtA.to_haplotypes()
    gtB = gt.take(pops[1], axis=1)
    htB = gtB.to_haplotypes()
    for hap1 in list(range(len(pops[0]))):
        for hap2 in list(combinations(range(len(pops[1])), 2)):
            ma = htA[:, [hap1]].count_alleles()
            mb = htB[:, hap2].count_alleles()
            jsfs = allel.joint_sfs(ma[:, 1], mb[:, 1])
            n1list.append(jsfs[0, 2] + jsfs[1, 0])
            n2list.append(jsfs[0, 1] + jsfs[1, 1])
            n3list.append(jsfs[0, 0] + jsfs[1, 2])
    c, k = estimCandK(n1list, n2list, n3list)
    return(c, k)


def countN1N2N3Fast(gt, pops):
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
            # fast approx
            k_hat = .75 * ((2*n3 + n2) / (n1 + n2 + n3))
            c_hat = (2*n3 - n2) / (2*n3 + n2)
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
                # fast approx
                k_hat = .75 * ((2*n3 + n2) / (n1 + n2 + n3))
                c_hat = (2*n3 - n2) / (2*n3 + n2)
                clist.append(c_hat)
                klist.append(k_hat)
    return(np.mean(clist), np.mean(klist))


def estimDiv(c, k, psmc, j):
    """Estimate divergence using eq 12
    """
    if psmc:
        try:
            assert j is True
        except AssertionError:
            print('PSMC require divTime option, defaulting to Constant')
            T_hat = -log(1-c)  # assumes constant popsize
        # check MSlatkins mathematica code for setting this up
        # check format of psmc, should be time in coalescent and then theta
        # msmc needs to be transformed to match psmc input
#        with open(psmc, 'r') as pw:
#            for line in pw:
#                pass
#        psmc_t = [0, 1, 2, 3, 4, 5]
#        psmc_N = [100, 200, 300, 500, 600]
#        prodNis = [exp(-((psmc_t[t+1] - psmc_t[t]) / 2*N)) for t, N in enumerate(psmc_N[:j])]
#        # not sure how to do numerical solving in python
#        T_hat = 1 - exp((T - psmc_t[j]) / 2*psmc_N[j]) * prodNis  # numerically solve for T
#        print("{}".format(T_hat/2*N0))
    else:
        T_hat = -log(1-c)  # assumes constant popsize
    return(T_hat)


def calcCI(gt, pops, psmc, j, boots):
    """
    """
    print("Running bootstraps...")
    T_hatlist = []
    # random resampling
    indices_rs = np.random.randint(0, len(gt), (1, boots, len(gt)))
    for b in range(boots):
        gt = gt.take(indices_rs[0][b], axis=0)
        c, k = countN1N2N3Fast(gt, pops)
        T_hatlist.append(estimDiv(c, k, psmc, j))
    # quantiles
    t_lowCI = np.percentile(T_hatlist, 0.025)
    t_highCI = np.percentile(T_hatlist, 0.975)
    return(t_lowCI, t_highCI)


if __name__ == "__main__":
    popset = args.pops
    psmc = args.piecewise
    j = args.divtime
    # pop_ix = [tuple(map(int, x)) for x in popset]
    pop_ix = [list(map(int, x)) for x in popset]
    if args.outgroup:
        outgroup_ix = list(map(int, args.outgroup))
    else:
        outgroup_ix = ''
    callset = loadvcf(args.vcfFile)
    gt = filterGT(callset, pop_ix, outgroup_ix)
    if args.mle:
        c, k = countN1N2N3MLE(gt, pop_ix)
    else:
        c, k = countN1N2N3Fast(gt, pop_ix)
    T_hat = estimDiv(c, k, psmc, j)
    t_LCI, t_HCI = calcCI(gt, pop_ix, psmc, j, args.boots)
    print("{} in 2Ne gens ({} - {})".format(T_hat, t_LCI, t_HCI))
