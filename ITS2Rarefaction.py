#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 10:30:37 2019

@author: scott
"""

from __future__ import print_function
from __future__ import division

import sys
import numpy as np


def rarefaction(Its2File):
    """samples with replacement from a blast file, depends on header
    """
    specieslist = []
    with open(Its2File) as its2:
        # header = its2.next()
        # spix = header.index("Species")
        for line in its2:
            # species = line.split()[spix]
            species = line.split()[5]
            specieslist.append(species)
    # resample and count
    seqs = range(10, len(species), 50)
    f = open("rarefaction.out", 'w')
    for s in seqs:
        rare = np.random.choice(specieslist, s)
        f.write("{}\t{}\n".format(s, len(set(rare))))
    f.close()
    return(None)


if __name__ == "__main__":
    rarefaction(sys.argv[1])
