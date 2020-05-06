#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 12:57:29 2019

@author: scott
"""

from __future__ import print_function
from __future__ import division
import argparse
from Bio import SeqIO
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--fasta", type=str, required=True, help="fasta file")
parser.add_argument('-anc', "--ancestral", type=str, help="position of the ancestral")
args = parser.parse_args()


def sfs(fastaFile, ancFile, step=1000):
    """
    """
    alignlist = []
    totalbp = 0
    thinbp = 0
    fasta_sequences = list(SeqIO.parse(fastaFile, 'fasta'))
    # load fasta
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
        alignlist.append(sequence)
    # load ancestral seq for polarization
    fasta_anc = list(SeqIO.parse(ancFile, 'fasta'))
    for fasta in fasta_anc:
        ancseq = str(fasta.seq)
    # calc sfs derived
    sfsarray = np.zeros(len(alignlist))
    sfsarray_thin = np.zeros(len(alignlist))
    thinStart = 0
    for bp in range(len(ancseq)):
        anc = ancseq[bp]
        al = [seq[bp] for seq in alignlist]
        if "N" in anc or "N" in al:
            continue
        else:
            alcnt = np.unique(al, return_counts=True)
            alleles = alcnt[0]
            freq = alcnt[1]
            if alleles.size > 2:
                continue
            elif alleles.size == 1:
                totalbp += 1
                thinbp += 1
            else:
                try:
                    ix = list(alcnt[0]).index(anc)
                    if ix == 0:
                        sfsarray[freq[1]] += 1
                        if bp - thinStart > step:
                            sfsarray_thin[freq[1]] += 1
                            thinStart = bp
                            thinbp += 1
                    else:
                        sfsarray[freq[0]] += 1
                        if bp - thinStart > step:
                            sfsarray_thin[freq[0]] += 1
                            thinStart = bp
                            thinbp += 1
                    totalbp += 1
                except ValueError:
                    continue
    print(totalbp, thinbp)
    return(sfsarray, sfsarray_thin)


if __name__ == "__main__":
    sfs, sfst = sfs(args.fasta, args.ancestral)
    print("{}".format(" ".join(map(str, sfs))))
    print("{}".format(" ".join(map(str, sfst))))
