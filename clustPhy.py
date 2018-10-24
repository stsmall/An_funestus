#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 17:28:42 2018
python clustPhy.py 100
@author: scott
"""

import glob
import sys

phyFile = glob.glob("*.phy")
clust = int(sys.argv[1])
lsorted = sorted(phyFile, key=lambda x: int(x.split(".")[3]))
i = 0
j = clust
while j < len(lsorted):
    clustFiles = lsorted[i:j]
    # get start
    start = open(clustFiles[0], 'r')
    start.readline()
    ls = start.readline().split()[0]
    s = ls.split("/")[-1].split("-")[0]
    start.close()
    # get end
    end = open(clustFiles[-1], 'r')
    end.readline()
    le = end.readline().split()[0]
    e = le.split("/")[-1].split("-")[1]
    end.close()
    # open file with proper coords
    f = open("bpp.{}-{}.txt".format(s, e), 'w')
    for phy in clustFiles:
        with open(phy, 'r') as p:
            header = p.next()
            f.write("\n")
            f.write(header)
            f.write("\n")
            for line in p:
                x = line.split()
                spname = x[0].split("-")[0]
                f.write("^{}{}{}\n".format(spname,' '*4, x[1]))
    f.close()
    i = j
    j += clust

# last window
j = len(lsorted)
clustFiles = lsorted[i:j]
# get start
start = open(clustFiles[0], 'r')
start.readline()
ls = start.readline().split()[0]
s = ls.split("/")[-1].split("-")[0]
start.close()
# get end
end = open(clustFiles[-1], 'r')
end.readline()
le = end.readline().split()[0]
e = le.split("/")[-1].split("-")[1]
end.close()
# open file with proper coords
f = open("bpp.{}-{}.txt".format(s, e), 'w')
for phy in clustFiles:
    with open(phy, 'r') as p:
        header = p.next()
        f.write("\n")
        f.write(header)
        f.write("\n")
        for line in p:
            x = line.split()
            spname = x[0].split("-")[0]
            f.write("^{}{}{}\n".format(spname,' '*4, x[1]))
f.write("\n")
f.close()
