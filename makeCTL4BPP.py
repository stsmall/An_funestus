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
    f = open("A01.bpp.{}.ctl".format(c), 'w')
    with open(sys.argv[1], 'r') as Sctl:
        for line in Sctl:
            if line.startswith("seqfile"):
                f.write("seqfile = {}\n".format(ctl))
            elif line.startswith("outfile"):
                f.write("outfile = {}.out\n".format(".".join(ctl.split(".")[0:-1])))
            elif line.startswith("mcmcfile"):
                f.write("mcmcfile = {}.msmc.txt\n".format(".".join(ctl.split(".")[0:-1])))
            else:
                f.write(line)
    c += 1
    f.close()
