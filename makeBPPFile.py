#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 16:23:36 2018
python makeBPPFile.py --gff gff --fasta FOO.fa [--exons] [--distance int] [--length int] [--chromlen int]
@author: scott
"""
from __future__ import print_function
from __future__ import division
import argparse
from Bio import SeqIO
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--gff", type=str, required=True, help="gff file")
parser.add_argument("--distance", type=int, default=2000, help="distance"
                    " between non-coding loci")
parser.add_argument("--length", type=int, default=1000, help="length for"
                    " non-coding loci")
parser.add_argument("--fasta", type=str, required=True, help="fasta file")
parser.add_argument("--clust", type=int, default=100, help="how many loci"
                    " to cluster into 1 file")
parser.add_argument("--prct", type=float, default=0.5)
parser.add_argument("--exons", action="store_true")
parser.add_argument("--chromlen", type=int, help="length of chromosome or contig")
args = parser.parse_args()


#TODO: test exons
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
            if not line.startswith("#"):
                x = line.split()
                chrom = x[0]
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
    return(cdsdict, chrom)


def getNonCDS(cdsdict, lengths, distance, exons, chromlen):
    """
    """
    noncdsdict = {}
    loci = 0
    for i, key in enumerate(cdsdict.keys()):
        k = "cds_" + str(i)
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
                    #TODO: CHROMLEN
                    next_start = chromlen
                    Sstart = end + distance
                    Send = Sstart + lengths
                    while next_start - Send > distance:
                        noncdsdict["ncds_" + str(loci)] = [int(Sstart), int(Send)]
                        Sstart = end + distance
                        Send = Sstart + lengths
                        loci += 1
                else:
                    break
        Sstart = end + distance
        Send = Sstart + lengths
        if next_start - Send > distance:
            noncdsdict["ncds_" + str(loci)] = [int(Sstart), int(Send)]
            loci += 1
    return(noncdsdict)


def bppFormatCDS(CDSdict, fastaFile, clust, exons, chrom, prct, just=10):
    """
    """
    fasta_sequences = list(SeqIO.parse(fastaFile, 'fasta'))
    # CDS
    skip_gaps = 0
    print("CDS file")
    loci = 0
    s = CDSdict["cds_0"][0]
    try:
        e = CDSdict["cds_{}".format(clust)][1]
    except KeyError:
        e = CDSdict["cds_{}".format(len(CDSdict.keys())-1)][1]
    out_file = open("CDS.bpp.{}.{}-{}.txt".format(chrom, s, e), 'w')
    for i in range(len(CDSdict.keys())):
        k = "cds_" + str(i)
        xlist = CDSdict[k]
        locuslist = []
        headerlist = []
        if loci > clust:
            out_file.close()
            try:
                s = CDSdict["cds_{}".format(i)][0]
                e = CDSdict["cds_{}".format(i+clust)][1]
            except KeyError:
                e = CDSdict["cds_{}".format(len(CDSdict.keys())-1)][1]
            out_file = open("CDS.bpp.{}.{}-{}.txt".format(chrom, s, e), 'w')
            loci = 0
        if exons:
            for fasta in fasta_sequences:
                headerlist.append(fasta.id)
                locus = ''
                for ex in xlist:
                    sequence = str(fasta.seq)
                    locus += sequence[xlist[0]:xlist[1]]
                locuslist.append(locus)
        else:
            firstREF = True
            gap = 0
            for fasta in fasta_sequences:
                header, sequence = fasta.id, str(fasta.seq)
                if firstREF:
                    if sequence[xlist[0]:xlist[1]].count("-") > 0:
                        while gap < sequence[xlist[0]:xlist[1]+gap].count("-"):
                            gap += sequence[xlist[0]:xlist[1]].count("-")
                locuslist.append(sequence[xlist[0]:xlist[1]+gap])
                headerlist.append(header)
                firstREF = False
        samples = len(headerlist)
        length = len(locuslist[0])
        # Ns check point
        try:
            if any(seqX.count("N")/length > prct for seqX in locuslist):
                # print("skipping, too many Ns")
                skip_gaps += 1
            else:
                out_file.write("\n{} {}\n\n".format(samples, length))
                for head, seq in zip(headerlist, locuslist):
                    out_file.write("^{}{}{}\n".format(head, ' '*(just-len(head)), seq))
                loci += 1
        except ZeroDivisionError:
            import ipdb;ipdb.set_trace()
    out_file.close()
    print("{} regions skipped due to excess gaps/N's".format(skip_gaps))
    return(None)


def bppFormatnCDS(nonCDSdict, fastaFile, clust, chrom, prct, just=10):
    """
    """
    skip_gaps = 0
    fasta_sequences = list(SeqIO.parse(fastaFile, 'fasta'))
    # nonCDS
    print("nonCDS file")
    loci = 0
    s = nonCDSdict["ncds_0"][0]
    try:
        e = nonCDSdict["ncds_{}".format(clust)][1]
    except KeyError:
        e = nonCDSdict["ncds_{}".format(len(nonCDSdict.keys())-1)][1]
    out_file = open("nCDS.bpp.{}.{}-{}.txt".format(chrom, s, e), 'w')
    for i in range(len(nonCDSdict.keys())):
        k = "ncds_" + str(i)
        xlist = nonCDSdict[k]
        locuslist = []
        headerlist = []
        if loci > clust:
            out_file.close()
            try:
                s = nonCDSdict["ncds_{}".format(i)][0]
                e = nonCDSdict["ncds_{}".format(i+clust)][1]
            except KeyError:
                s = nonCDSdict["ncds_{}".format(i)][0]
                e = nonCDSdict["ncds_{}".format(len(nonCDSdict.keys())-1)][1]
            out_file = open("nCDS.bpp.{}.{}-{}.txt".format(chrom, s, e), 'w')
            loci = 0
        for fasta in fasta_sequences:
            header, sequence = fasta.id, str(fasta.seq)
            locuslist.append(sequence[xlist[0]:xlist[1]])
            headerlist.append(header)
        samples = len(headerlist)
        length = len(locuslist[0])
        # Ns check point
        if any(seqX.count("N")/length > prct for seqX in locuslist) or any(seqX.count("-")/length > prct for seqX in locuslist):
            # print("skipping, too many Ns/gaps")
            skip_gaps += 1
        else:
            out_file.write("\n{} {}\n\n".format(samples, length))
            for head, seq in zip(headerlist, locuslist):
                out_file.write("^{}{}{}\n".format(head, ' '*(just-len(head)), seq))
            loci += 1
    out_file.close()
    print("{} regions skipped due to excess gaps/N's".format(skip_gaps))
    return(None)


if __name__ == "__main__":
    gffFile = args.gff
    length = args.length
    exons = args.exons
    prctmiss = args.prct
    assert exons is False
    distance = args.distance
    fastaFile = args.fasta
    clust = args.clust
    CDSdict, chrom = getCDS(gffFile, exons)
    nonCDSdict = getNonCDS(CDSdict, length, distance, exons, args.chromlen)
    bppFormatCDS(CDSdict, fastaFile, clust, exons, chrom, prctmiss, just=10)
    bppFormatnCDS(nonCDSdict, fastaFile, clust, chrom, prctmiss, just=10)
