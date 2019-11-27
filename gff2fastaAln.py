#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 12:21:45 2019
@author: Scott T. Small

This module creates input files for the methods of Thawornattana 2018.
Using a gff file it parses out individual alignments of exons and non-coding
sequences.

Example
-------

    $ python gff2fastaAln.py --gff FOO.gff --aln FOO.aln.fa [--exons] [--distance int] [--length int] [--chromlen int]

Notes
-----

    This is an example of an indented section. It's like any other section,
    but the body is indented to help it stand out from surrounding text.

If a section is indented, then a section break is created by
resuming unindented text.

"""
from typing import List, Set, Dict, Tuple, Optional
from Bio import SeqIO
from dataclasses import dataclass
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--gff", type=str, help="gff file")
parser.add_argument("--aln", type=str, help="alignment in fasta"
                    " format")
parser.add_argument("--minCDS", type=int, default=100, help="min lenght for including"
                    " coding loci")
parser.add_argument("--prct", type=float, default=0.5, help="percent sequence"
                    " length to allow that is Ns or gaps")
parser.add_argument("--distance", type=int, default=2000, help="min distance"
                    " between non-coding loci")
parser.add_argument("--mxlength", type=int, default=1000, help="max length for"
                    " non-coding loci")
parser.add_argument("--mnlength", type=int, default=100, help="min length for"
                    " non-coding loci")
parser.add_argument("--chromlen", type=int, help="length of chromosome or contig")
parser.add_argument("--bpp", action="store_true", help="bpp input files")
parser.add_argument("--clust", type=int, default=100, help="how many loci"
                    " to cluster into 1 file")
args = parser.parse_args()


@dataclass
class CDS:
    start: int
    end: int


@dataclass
class nCDS:
    start: int
    end: int


def getCDS(gffFile, minCDS, key="CDS"):
    """Converts gff file to dict with coordinates of each CDS

    Parameters
    ----------
    gffFile: File
        file in gff or gff3 format
    minCDS: int
        min length of coding loci

    Returns
    -------
    cds: dict
        dict with values of dataclass
    chrom: str
        value of chromosome

    """
    cds = {}  # type: Dict
    i = 0  # type: int
    with open(gffFile, 'rt') as gff:
        for line in gff:
            if not line.startswith("#"):
                x = line.split()
                chrom = x[0]  # type: str
                feature = x[2]  # type: str
                if key in feature:
                    start = int(x[3])  # type: int
                    end = int(x[4])  # type: int
                    if start < end:
                        start = int(x[4])
                        end = int(x[3])
                    if end - start >= minCDS:
                        cds["cds_" + str(i)] = CDS(start-1, end)
                        i += 1
    return(cds, chrom)


def getNonCDS(cdsdict, mxlen: int, mnlen: int, distance: int, chromlen: int):
    """Returns the non-coding loci around the coding sequence

    Parameters
    ----------
    cdsdict: dict
        dict with values of dataclass: CDS.start, CDS.stop
    mxlen: int
        maximum length of non-coding loci
    mnlen: int
        min length of non-coding loci
    distance: int
        distance between non-coding loci
    chromlen: int
        length of chromosome

    Returns
    -------
    non_cds: dict
        dict with values of dataclass

    """
    non_cds = {}  # type: Dict
    loci = 0  # type: int
    for i in range(0, len(cdsdict)):
        start = cdsdict["cds_" + str(i)].end  # type: int
        end = cdsdict["cds_" + str(i+1)].start  # type: int
        nstart = start
        nend = nstart + mxlen
        if nend < end:
            while nend < end:
                non_cds["ncds_"+str(loci)] = nCDS(nstart, nend)
                nstart += distance
                nend = nstart + mxlen
                loci += 1
            if (end - nstart) >= mnlen:
                non_cds["ncds_"+str(loci)] = nCDS(nstart, end)
                loci += 1
        else:
            if (end - start) >= mnlen:
                non_cds["ncds_"+str(loci)] = nCDS(start, end)
                loci += 1

    # end of chrom
    nstart = end
    nend = end + mxlen
    if nend < chromlen:
        while nend < chromlen:
            non_cds["ncds_"+str(loci)] = nCDS(nstart, nend)
            nstart += distance
            nend = nstart + mxlen
            loci += 1
        if (chromlen - nstart) >= mnlen:
            non_cds["ncds_"+str(loci)] = nCDS(nstart, end)
    else:
        if (chromlen - end) >= mnlen:
            non_cds["ncds_"+str(loci)] = nCDS(end, chromlen)
    return(non_cds)


def formatFasta(fname, gff_dict, fastaFile, clust: int, chrom: str, prct: float, bpp: bool, just=10):
    """Creates BPP input files or if clust=1 and BPP=False, then a fasta file
    for each CDS and non-CDS.

    Parameters
    ----------
    fname: str
        cds or ncds
    gff_dict: dict
        dict with values of dataclass: fname.start, fname.stop
    fastaFile: file
        alignment in fasta format
    clust: int
        number of loci to push into 1 file
    chrom: str
        chromosome name
    prct: float
        cutoff for percent of missing data in alignment
    bpp: bool
        format for BPP True or False
    just: int
        justified for BPP

    Returns
    -------
    None

    """

    fasta_sequences = list(SeqIO.parse(fastaFile, 'fasta'))
    skip_gaps = 0
    f"{fname} file"
    loci = 0
    s = gff_dict[f"{fname}_0"].start
    try:
        e = gff_dict[f"{fname}_{clust-1}"].end
    except KeyError:
        e = gff_dict[f"{fname}_{len(gff_dict)-1}"].end
    out_file = open(f"CDS.bpp.{chrom}.{s}-{e}.txt", 'w')
    for i in range(len(gff_dict)):
        k = [f"{fname}_{str(i)}"]
        locuslist = []
        headerlist = []
        if loci >= clust:
            out_file.close()
            try:
                s = gff_dict[k].start
                e = gff_dict[f"{fname}_{i + clust-1)}"].end
            except KeyError:
                e = gff_dict[f"{fname}_{len(gff_dict)-1)}"].end
            out_file = open(f"{fname}.bpp.{chrom}.{s}-{e}.txt", 'w')
            loci = 0
        else:
            for fasta in fasta_sequences:
                header, sequence = fasta.id, str(fasta.seq)
                locuslist.append(sequence[gff_dict[k].start:gff_dict[k].end])
                headerlist.append(header)
        samples = len(headerlist)
        seqlen = len(locuslist[0])
        # Ns check point
        try:
            if any((seqX.count("N")/seqlen) > prct for seqX in locuslist):
                # print("skipping, too many Ns")
                skip_gaps += 1
            else:
                out_file.write(f"{samples} {seqlen}\n\n")
                for head, seq in zip(headerlist, locuslist):
                    if bpp:
                        out_file.write(f"^{head}{' '*(just-len(head))}{seq}\n")
                    else:
                        out_file.write(f">{head}\n{seq}\n")
                loci += 1
        except ZeroDivisionError:
            import ipdb;ipdb.set_trace()
    out_file.close()
    print(f"{skip_gaps} regions skipped due to excess gaps/N's")
    return(None)


if __name__ == "__main__":
    gffFile = args.gff
    length = args.length
    prctmiss = args.prct
    distance = args.distance
    fastaFile = args.fasta
    clust = args.clust
    cds_dict, chrom = getCDS(gffFile)
    ncds_dict = getNonCDS(cds_dict, length, distance, args.chromlen)
    formatFasta("cds", cds_dict, fastaFile, clust, chrom, prct, args.bpp, just=10)
    formatFasta("ncds", ncds_dict, fastaFile, clust, chrom, prct, args.bpp, just=10)
