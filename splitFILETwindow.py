#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 16:02:46 2019

@author: scott
"""

import sys

f = open("{}.sp1.bed".format(sys.argv[1]), 'w')
g = open("{}.sp2.bed".format(sys.argv[1]), 'w')
with open(sys.argv[1], 'r') as sig:
    for line in sig:
        x = line.split()
        sp1_count = x[3].count("1")
        sp2_count = x[3].count("2")
        if sp1_count > sp2_count:
            f.write(line)
        else:
            g.write(line)
