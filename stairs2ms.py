#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 16:42:04 2018
stairway2ms.py -i stairway.final.summary -L 1000
@author: scott
"""

from __future__ import print_function
from __future__ import division
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i', "--inFile", type=str, required=True, help="stairway file")
parser.add_argument('-L', "--LocusLen", type=int, default=1000, help="locus length")
parser.add_argument('-r', "--rhorat", type=float, help="ratio of mu/rho")
parser.add_argument('-e', "--epochs", type=int, help="epochs to print")
parser.add_argument("--subpop",type=int, default=0)
parser.add_argument("--time",type=int, default=0)
args = parser.parse_args()


def stairs2ms(inFile, locLen, rhorat, epochs, subpop, time):
    """
    """
    msmc2ms2 = []
    with open(inFile,  'r') as stair:
        for line in stair:
            if line.startswith("mutation"):
                line = next(stair)
                if time > 0:
                    gens = float(line.split()[5])
                    while time > gens:
                        line = next(stair)
                        gens = float(line.split()[5])
                x = line.split()
                theta0M = float(x[2])
                theta0L = float(x[3])
                theta0U = float(x[4])
                Ne0 = float(x[6])
                for line in stair:
                    x = line.split()
                    theta = float(x[2])
                    thetaL = float(x[3])
                    thetaU = float(x[4])
                    gens = float(x[5])
                    Ne = float(x[6])
                    msmc2ms2.append([gens/(4*Ne0), subpop, theta/theta0M])
    if epochs:
        reducedEpochs = int(len(msmc2ms2)/epochs)
        msmc2ms = msmc2ms2[0::reducedEpochs]
    else:
        msmc2ms = msmc2ms2
    tM = theta0M * locLen
    tL = theta0L * locLen
    tU = theta0U * locLen
    size = ' -en '.join(' '.join(str(x) for x in timestep) for timestep in msmc2ms)
    print("-Pt {} {}".format(tL, tU))
    print("-Pr {} {}".format(tL*rhorat, tU*rhorat))
    print("-t {} -r {} -en {}".format(tM, tM * rhorat, size))
    return(None)


if __name__ == "__main__":
    stairs2ms(args.inFile, args.LocusLen, args.rhorat, args.epochs, args.subpop, args.time)
