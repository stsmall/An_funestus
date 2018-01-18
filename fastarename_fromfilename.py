#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 13:19:27 2017

@author: scott
"""
import glob

fasta_files = glob.glob("*.fasta")

for fa in fasta_files:
    with open(fa + "_rename", 'w') as out:
        with open(fa, 'r') as fa_in:
            for line in fa_in:
                if line.startswith(">"):
                    out.write(">{}\n".format("_".join(fa.split(".")[0:3])))
                else:
                    out.write(line)
