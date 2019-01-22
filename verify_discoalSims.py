#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 14:46:17 2019

@author: scott
"""
import glob

files = glob.glob("*msOut")
for ms in files:
    f = open("{}2".format(ms), 'w')
    with open(ms, 'r') as sims:
        line = sims.next()
        samples = int(line.split()[1])
        f.write(line)
        for line in sims:
            if line.startswith("discoal"):
                continue
            else:
                if line.startswith("//"):
                    segline = sims.next()
                    posline = sims.next()
                    scount = 0
                    hapline = []
                    while line.strip():
                        hapline.append(line)
                        scount += 1
                        line = sims.next()
                    if scount == samples:
                        f.write("//\n")
                        f.write(segline)
                        f.write(posline)
                        for h in hapline:
                            f.write(h)
                        f.write("\n")
                else:
                    f.write(line)  # this should be blank lines
    f.close()
