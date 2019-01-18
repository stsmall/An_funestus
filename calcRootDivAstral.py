#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 16:48:42 2019

@author: scott
"""

import re
import numpy as np
import sys

x = []
with open(sys.argv[1], 'r') as tr:
    for line in tr:
        t = re.findall(r'\d?\.?\d+\.?\d?(?=:)|(?<=\):)\d+\.?\d+', line.strip())
        x.append(sum(map(float, t)))
        if sum(map(float, t)) > 10:
            print(line)
print(np.mean(x))
print(np.quantile(x, 0.025))
print(np.quantile(x, 0.975))
