#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 13:39:32 2018

@author: scott
"""

import sys
import gzip

fastaFile = sys.argv[1]
with gzip.open("{}.1.fa.gz".format(fastaFile.rstrip(".fa.gz")), 'wb') as f:
    if "gz" in fastaFile:
        with gzip.open(fastaFile, 'rb') as fasta:
            line = next(fasta)
            f.write(line)
            for line in fasta:
                if line.startswith(b">"):
                    break
                else:
                    f.write(line)
        f.close()
    else:
        with open(fastaFile, 'r') as fasta:
            line = next(fasta)
            f.write(line)
            for line in fasta:
                if line.startswith(">"):
                    break
                else:
                    f.write(line)
        f.close()
