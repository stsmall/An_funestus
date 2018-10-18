#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 16:23:36 2018
python makeBPPfile.py --gff gff.bed --distance 2000 --length 1000 --fasta FOO.fa
@author: scott
"""

import argparse
from Bio import SeqIO
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--gff", type=str, required=True, help="bed file from gff")
parser.add_argument("--distance", type=int, default=2000, help="distance"
                    " between non-coding loci")
parser.add_argument("--length", type=int, default=1000, help="length for"
                    " non-coding loci")
parser.add_argument("--fasta", type=str, required=True, help="fasta file")
parser.add_argument("--clust", type=int, default=1000000, help="how many loci"
                    " to cluster into 1 file")
parser.add_argument("--exons", action="store_true")
parser.add_argument("--chromlen", type=int, help="length for non-coding loci")
args = parser.parse_args()


def getCDS(gffFile, exons):
    """
    """
    if exons:
        cdsdict = defaultdict(list)
    else:
        cdsdict = {}
    i = 0
    with open(gffFile, 'r') as gff:
        for line in gff:
            x = line.split()
            feature = x[2]
            start = x[3]
            end = x[4]
            if exons:
                if "CDS" in feature:
                    while "CDS" in feature:
                        cdsdict["cds_" + str(i)].append((int(start)-1, int(end)))
                        line = gff.next()
                        x = line.split()
                        feature = x[2]
                        start = x[3]
                        end = x[4]
                    i += 1
            else:
                if "gene" in feature:
                    cdsdict["cds_" + str(i)] = [int(start)-1, int(end)]
                    i += 1
    return(cdsdict)


def getNonCDS(cdsdict, lengths, distance, exons, chromlen):
    """
    """
    noncdsdict = {}
    for i, k in enumerate(cdsdict.keys()):
        if i == 0:
            end = 0
            if exons:
                next_start = cdsdict[k][0][0]
            else:
                next_start = cdsdict[k][0]
        else:
            try:
                if exons:
                    end = cdsdict[k][-1][1]
                    next_start = cdsdict["cds_" + str(i+1)][0][0]
                else:
                    end = cdsdict[k][1]
                    next_start = cdsdict["cds_" + str(i+1)][0]
            except KeyError:
                # beyond the genes in gff
                if chromlen:
                    next_start = chromlen
                else:
                    break
        Sstart = end + distance
        Send = Sstart + lengths
        while next_start - Send > distance:
            noncdsdict["ncds" + str(i)] = [int(Sstart)-1, int(Send)]
            Sstart = end + distance
            Send = Sstart + lengths
    return(noncdsdict)


def bppFormat(CDSdict, nonCDSdict, fastaFile, clust, exons):
    """
    """
    fasta_sequences = SeqIO.parse(fastaFile, 'fasta')
    # CDS
    loci = 1
    out_file = open("CDS.bpp.{}.txt".format(loci), 'w')
    for cds in CDSdict:
        xlist = CDSdict[cds]
        locuslist = []
        headerlist = []
        if loci > clust:
            out_file.close()
            out_file = open("CDS.bpp.{}.txt".format(loci), 'w')
        if exons:
            for fasta in fasta_sequences:
                headerlist.append(fasta.id)
                locus = ''
                for ex in xlist:
                    sequence = str(fasta.seq)
                    locus += sequence[xlist[0]:xlist[1]]
                locuslist.append(locus)
        else:
            for fasta in fasta_sequences:
                header, sequence = fasta.id, str(fasta.seq)
                locuslist.append(sequence[xlist[0]:xlist[1]])
                headerlist.append(header)
        samples = len(headerlist)
        length = len(locuslist[0])
        out_file.write("\n{} {}\n\n".format(samples, length))
        for head, seq in zip(headerlist, locuslist):
            out_file.write("^{} {}\n".format(head, seq))
        loci += 1
    out_file.close()

    # nonCDS
    loci = 1
    out_file = open("nCDS.bpp.{}.txt".format(loci), 'w')
    for cds in nonCDSdict:
        xlist = nonCDSdict[cds]
        locuslist = []
        headerlist = []
        if loci > clust:
            out_file.close()
            out_file = open("nCDS.bpp.{}.txt".format(loci), 'w')
        for fasta in fasta_sequences:
            header, sequence = fasta.id, str(fasta.seq)
            locuslist.append(sequence[xlist[0]:xlist[1]])
            headerlist.append(header)
        samples = len(headerlist)
        length = len(locuslist[0])
        out_file.write("\n{} {}\n\n".format(samples, length))
        for head, seq in zip(headerlist, locuslist):
            out_file.write("^{} {}\n".format(head, seq))
        loci += 1
    out_file.close()
    return(None)


if __name__ == "__main__":
    gffFile = args.gff
    length = args.length
    exons = args.exons
    distance = args.distance
    fastaFile = args.fasta
    clust = args.clust
    CDSdict = getCDS(gffFile, exons)
    nonCDSdict = getNonCDS(CDSdict, length, distance, exons, args.chromlen)
    bppFormat(CDSdict, nonCDSdict, fastaFile, clust, exons)
