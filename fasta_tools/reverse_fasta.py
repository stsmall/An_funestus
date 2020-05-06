#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 17:11:56 2017

@author: scott
"""

from pyfaidx import Fasta
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--fasta", type=str,
                    help='fasta file')
parser.add_argument('-c', "--coordinates", type=str,
                    help='coordinates file')
args = parser.parse_args()


def loadfasta(fastafile):
    """
    """
    ff = Fasta(fastafile, mutable=True)
    return(ff)


def loadcoords(coords):
    """
    """
    coorddict = {}
    with open(coords, 'r') as c:
        for line in c:
            x = line.strip().split()
            coorddict[x[2]] = (int(x[0] - 1), int(x[1]))
    return(coorddict)


def reversefasta(coord, fasta):
    """
    """
    for contig in coord.keys():
        s = coord[contig][0]
        e = coord[contig][1]
        rev = fasta[contig][s:e].reverse.seq
        fasta[contig][s:e] = rev
    return(None)


if __name__ == "__main__":
    fastafile = args.fasta
    coord = args.coordinates
    ff = loadfasta(fastafile)
    cdict = loadcoords(coord)
    reversefasta(cdict, ff)

    # CHR3
    fastafile = "Anfunestus.3.fa"
    phase = Fasta(fastafile, mutable=True)
    x = phase["PGA_scaffold1_scaffold2__421_contigs__length_66724930_3L_3R"][46586560:66001209].reverse.complement.seq
    phase["PGA_scaffold1_scaffold2__421_contigs__length_66724930_3L_3R"][46586560:66001209] = x

    # CHR2
    fastafile = "Anfunestus.2.fa"
    phase = Fasta(fastafile, mutable=True)
    # 2L
    f1 = phase["PGA_scaffold0__405_contigs__length_98997778_2L_2R"][70863710:74366026].reverse.complement.seq
    phase["PGA_scaffold0__405_contigs__length_98997778_2L_2R"][70863710:74366026] = f1
    f2 = phase["PGA_scaffold0__405_contigs__length_98997778_2L_2R"][72881800:98204790].reverse.complement.seq
    phase["PGA_scaffold0__405_contigs__length_98997778_2L_2R"][72881800:98204790] = f2

    f3 = phase["PGA_scaffold0__405_contigs__length_98997778_2L_2R"][54617112:70863771].reverse.complement.seq
    phase["PGA_scaffold0__405_contigs__length_98997778_2L_2R"][54617112:70863771] = f3

   #3R = phase["3L_3R"][46586561:].reverse
   #3La = phase["3L_3R"][0:46586455]
   #2R = phase["2L_2R"][0:54617112]
   #2L = phase["2L_2R"][54617112:].reverse