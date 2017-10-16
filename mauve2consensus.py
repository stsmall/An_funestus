#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function
"""
Created on Mon Oct 16 13:12:38 2017
@author: stsmall

I made this script to create a consensus between 2 de novo genome assemblies
each from a difference individual but from the same species. The de novo
assemblies were then scaffolded in ragout. N's in 1 scaffolded assembly were
somtimes present as sequence in the other. If you align the 2 genomes with
Mauve to product a xmfa, you can then produce a consensus with fewer gaps than
either alone
"""
from Bio import AlignIO
from pyfaidx import Fasta
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', "--fasta", type=str, required=True,
                    help='fasta for base')
parser.add_argument('-r1', "--ref1", type=str, required=True,
                    help='aligned ref 1')
parser.add_argument('-r2', "--ref2", type=str, required=True,
                    help='aligned ref 2')
parser.add_argument('-r', "--ref", type=str, required=True,
                    help='aligned ref')
parser.add_argument('-x', "--xmfa", type=str, required=True,
                    help='mauve xmfa')
args = parser.parse_args()


def makesense(alnarr, header, r1, r2, ref):
    """
    """
    print("making sense...")
    gap = 0
    gapfill = 0
    gapfillmask = 0
    sense = ""
    h = [hd.split(".")[0] for hd in header]
    seq = h.index(r1)
    seq1 = alnarr[seq]
    seq1 = np.char.upper(seq1)
    if len(header) == 3:
        alt = h.index(r2)
        seq2 = alnarr[alt]
        seq2 = np.char.upper(seq2)
        seq3 = alnarr[h.index(ref)]
        seq3 = np.char.lower(seq3)
        for i, base in enumerate(seq1):
            base1 = base
            base2 = seq2[i]
            base3 = seq3[i]
            if base1 == "N":
                gap += 1
                if base2 != "-" and base2 != "N":
                    sense += base2  # fill N with base from other seq
                    gapfill += 1
                elif base3 != "-" and base2 != "N":
                    sense += base3
                    gapfillmask += 1
                else:
                    sense += base1
            elif base1 != "-":
                sense += base1  # keep other bases
            else:
                pass
    if len(header) == 2:
        if r2 in h:
            alt = h.index(r2)
            seq2 = alnarr[alt]
            seq2 = np.char.upper(seq2)
        else:
            seq2 = alnarr[h.index(ref)]
            seq2 = np.char.lower(seq2)
        for i, base in enumerate(seq1):
            base1 = base
            base2 = seq2[i]
            if base1 == "N":
                gap += 1
                if base2 != "-" and base2 != "N":
                    sense += base2  # fill N with base from other seq
                    gapfill += 1
                else:
                    sense += base1
            elif base1 != "-":
                sense += base1  # keep other bases
            else:
                pass  # no point in keeping alignment gaps
    fill = gap - gapfill
    print("gaps filled: {}, {}".format(fill, gapfillmask))
    return(sense, fill)


def parsexmfa(xmfa, r1, r2, ref):
    """
    """
    print("parsing...")
    r1 = r1
    r2 = r2
    consensusdict = {}
    gapfill = 0
    alignment = AlignIO.parse(open(xmfa), "mauve")
    for aln in alignment:  # each alignment block
        header = []
        if len(aln) > 1:
            for record in aln:
                header.append(record.id)
            for rec in header:
                if r1 in rec:
                    pos = rec.split("/")[1]
                    pos = pos.replace("-", ":")
                    alignarr = np.array([list(r) for r in aln], np.character)
                    sense, fill = makesense(alignarr, header, r1, r2, ref)
                    gapfill += fill
                    consensusdict[pos] = sense
    print("total gaps filled:{}".format(gapfill))
    return(consensusdict)


def fillgaps(consensusdict, fasta):
    """
    """
    print("filling consensus...")
    fastascaf = Fasta(fasta, mutable=True)
    for chrom in fastascaf.keys():
        for s in consensusdict.keys():
            t1 = int(s.split(":")[0])
            t2 = int(s.split(":")[1])
            assert (t2 - t1) == len(fastascaf[chrom][t1:t2].seq)
            fastascaf[chrom][t1:t2] = consensusdict[s]
    return(None)


if __name__ == "__main__":
    fasta = args.fasta
    r1 = args.ref1
    r2 = args.ref2
    ref = args.ref
    xmfa = args.xmfa
    condict = parsexmfa(xmfa, r1, r2, ref)
    fillgaps(condict, fasta)
