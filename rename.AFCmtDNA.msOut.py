# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import re
import sys

d = {"1": "Fun1",
     "2": "Fun2",
     "3": "Fun3",
     "4": "Fun4",
     "5": "Fun5",
     "6": "Fun6",
     "7": "Fun7",
     "8": "Fun8",
     "9": "Fun9",
     "10": "Fun10",
     "11": "Fun11",
     "12": "Fun12",
     "13": "Fun13",
     "14": "Fun14",
     "15": "Fun15",
     "16": "Lik1",
     "17": "Lik2",
     "18": "Lik3",
     "19": "Lik4",
     "20": "Lik5",
     "21": "Van1",
     "22": "Van2",
     "23": "Van3",
     "24": "Van4",
     "25": "Van5",
     "26": "Van6",
     "27": "Van7",
     "28": "Van8",
     "29": "Van9",
     "30": "Van10",
     "31": "Lon1",
     "32": "Lon2",
     "33": "Lon3",
     "34": "Lon4",
     "35": "Lon5",
     "36": "Lon6",
     "37": "Lon7",
     "38": "Lon8",
     "39": "Par1",
     "40": "Par2",
     "41": "Par3",
     "42": "Par4",
     "43": "Par5",
     "44": "Par6",
     "45": "Par7",
     "46": "Par8",
     "47": "Par9",
     "48": "Par10"}

f = open("mtDNA.renamed.out", 'w')
with open(sys.argv[1], 'r') as tre:
    for line in tre:
        for k in d.keys():
            p1 = re.compile(fr",({k}):")
            p2 = re.compile(fr"\(({k}):")
            m1 = re.findall(p1, tre)
            m2 = re.findall(p2, tre)
            if len(m1) == 1:
                tre = re.sub(p1, f",{d[k]}:", tre)
            elif len(m2) == 1:
                tre = re.sub(p2, f"({d[k]}:", tre)
        f.write(tre)
