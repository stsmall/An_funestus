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
from itertools import product
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, required=True,
                    help="vcf file of variants")
parser.add_argument('-p', "--pops", nargs='+', action="append",
                    required=True, help="index of populations"
                    "-p 1 2 3 -g 6 7 8 -p 11 12 13")
parser.add_argument('-o', "--outgroup", required=True,
                    help="index of outgroup")
parser.add_arguement('-psmc', "--piecewise", help="use psmc for population"
                     "size changes, if None then assumes constant")
parser.add_arguement('-N0', "--effectivesize", help="pop effective size")
parser.add_arguement('-j', "--divtime", help="a boundary where j < Tdiv < j+1"
                     "in coalescent units")
parser.add_arguement("--mle", action="store_true", help="use MLE, default is"
                     "fast approx. MLE requires larger popsizes")
args = parser.parse_args()


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


def countPatternsMLE(vcfFile, pops, outgroup):
    """Count patterns from VCF
    """
    # use allel to load a vcf, then count patterns of n1,n2,n3 for each
    # combination and estimating k1 and c.
    for combinations in pops:  # 1 from popA, 2 from popB
        # n1 S/SS, s/ss
        n1list.append(n1)
        # n2 S/Ss, s/Ss
        n2list.append(n2)
        # n3 S/ss, s/SS
        n3list.append(n3)
    c, k = estimCandK(n1, n2, n3)
    return(c, k)


def countPatternsFast(vcfFile, pops, outgroup):
    """Count patterns from VCF
    """
    # use allel to load a vcf, then count patterns of n1,n2,n3 for each
    # combination and estimating k1 and c.
    clist = []
    klist = []
    for combinations in pops:  # 1 from popA, 2 from popB
        # n1 S/SS, s/ss
        n1 =
        # n2 S/Ss, s/Ss
        n2 =
        # n3 S/ss, s/SS
        n3 =
        # fast approx
        k_hat = .75 * ((2*n3m + n2m) / (n1m + n2m + n3m))
        c_hat = (2*n3m - n2m) / (2*n3m + n2m)
        clist.append(c_hat)
        klist.append(k_hat)
    return(np.mean(clist), np.mean(klist))


def estimDiv(c, k, psmc, j, N0):
    """Estimate divergence using eq 12
    """
    if psmc:
        # check format of psmc, should be time in coalescent and then theta
        # msmc needs to be transformed to match psmc input
        with open(psmc, 'r') as pw:
            for line in pw:
                pass
        psmc_t = [0, 1, 2, 3, 4, 5]
        psmc_N = [100, 200, 300, 500, 600]
        prodNis = [exp(-((psmc_t[t+1] - psmc_t[t]) / 2*N)) for t, N in enumerate(psmc_N[:j])]
        T_hat = 1 - exp((T - psmc_t[j]) / 2*psmc_N[j]) * prodNis  # numerically solve for T
        print("{}".format(T_hat/2*N0))
    else:
        T_hat = -2*N0*log(1-c)  # assumes constant popsize
        print("{}".format(T_hat))
    return(T_hat)


if __name__ == "__main__":
    popset = args.pops
    pop_ix = [list(map(int, x)) for x in popset]
    outgroup_ix = map(int, args.outgroup)
    if args.mle:
        c, k = countPatternsMLE(args.vcfFile, pop_ix, outgroup_ix)
    else:
        c, k = countPatternsFast(args.vcfFile, pop_ix, outgroup_ix)
    T_hat = estimDiv(c, k, args.piecewise, args.divtime, args.effectivesize)
