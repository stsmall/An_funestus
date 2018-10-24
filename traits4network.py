#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 19:04:26 2018

@author: scott
"""

import sys
header = ["barbirostris", "nitidus", "peditaeniatus", "unknown"]
f = open("t", 'w')
with open(sys.argv[1], 'r') as log:
    log.next()  # skip header
    for line in log:
        try:
            sp = line.split()[0]
            headerlist = [0, 0, 0, 0]
            headerlist[header.index(sp.split("_")[0])] += 1
            line = log.next()
            while line.startswith("\t"):
                spix = header.index(line.split()[0].split("_")[0])
                headerlist[spix] += 1
                line = log.next()
            f.write("{},{}\n".format(sp, ",".join(map(str, headerlist))))
        except StopIteration:
            break
f.close()
