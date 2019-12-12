#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 17:34:46 2018

@author: scott
"""
import sys
import glob
ctl_files = glob.glob("*.txt")
c = 1
for ctl in ctl_files:
    coord = ctl.split(".")[2]
    chrom = ctl.split("."[1])
    cds = ctl.split(".")[0]
    counts = ctl.split(".")[3]
    f = open(f"A01.bpp.{cds}.{chrom}.{coord}.ctl", 'w')
    with open(sys.argv[1], 'r') as Sctl:
        for line in Sctl:
            if line.startswith("seqfile"):
                f.write("seqfile = {}\n".format(ctl))
            elif line.startswith("outfile"):
                f.write("outfile = {}.out\n".format(".".join(ctl.split(".")[0:-1])))
            elif line.startswith("nloci"):
                f.write(f"nloci = {counts}\n")
            elif line.startswith("mcmcfile"):
                f.write("mcmcfile = {}.mcmc.txt\n".format(".".join(ctl.split(".")[0:-1])))
            else:
                f.write(line)
    c += 1
    f.close()
