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
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', "--fasta", type=str, required=True,
                    help='fasta for base')
parser.add_argument('-r1', "--ref1", type=str, required=True,
                    help='aligned ref 1')
parser.add_argument('-r2', "--ref2", type=str, required=True,
                    help='aligned ref 2')
parser.add_argument('-x', "--xmfa", type=str, required=True,
                    help='mauve xmfa')
args = parser.parse_args()


def makesense(seq1, seq2):
    """
    """
    gap = 0
    gapfill = 0
    sense = ""
    for i, base in enumerate(seq1):
        base1 = base.upper()
        base2 = seq2[i].upper()
        if base1 == base2:
            sense += base1
        else:
            if base1 == "N" and base2 == "N":
                sense += base1
            elif base1 == "-":
                gap += 1
                if base2 == "-":
                    sense += base1
                else:
                    sense += base2
                    gapfill += 1
            else:
                sense += base1
    print("gaps filled: {}".format(gap - gapfill))
    return(sense)


def parsexmfa(xmfa, r1, r2):
    """
    """
    r1 = r1
    r2 = r2
    consensusdict = {}
    alignment = AlignIO.parse(open(xmfa), "mauve")
    for aln in alignment:  # each alignment block
        if len(aln) > 1:
            for record in aln:
                header = record.id
                if r1 in header:
                    pos = record.id.split("/")[1]
                    pos = pos.replace("-", ":")
                    seq1 = record.seq
                elif r2 in header:
                    seq2 = record.seq
                else:
                    pass
            sense = makesense(seq1, seq2)
            consensusdict[pos] = sense
    return(consensusdict)


def fillgaps(consensusdict, fasta):
    """
    """
    fastascaf = Fasta(fasta, mutable=True)
    for chrom in fastascaf.keys():
        for s in consensusdict.keys():
            t1 = int(s.split(":")[0])
            t2 = t1 + 10
            print(fastascaf[chrom][t1:t2].seq)
            fastascaf[chrom][int(s)] = consensusdict[s]
            print(fastascaf[chrom][t1:t2].seq)
    return(None)


if __name__ == "__main__":
    fasta = args.fasta
    r1 = args.ref1
    r2 = args.ref2
    xmfa = args.xmfa
