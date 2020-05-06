#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 12:47:16 2018

@author: scott
"""
from Bio import SeqIO
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', "--fasta", required=True, type=str,
                    help='path to fasta file')
parser.add_argument('-aa', '--ancestral', required=True, type=str,
                    help='ancestral fasta')
args = parser.parse_args()


def vcf2fasta(fastaFile, ancAllele):
    """reads in fasta, changes base using dictionary entry. note that fasta
       will be 0 based so position wll be -1
       notes on Biopython: SeqIO.to_dict() which builds all sequences into a
       dictionary and save it in memory
       SeqIO.index() which builds a dictionary without putting the sequences
       in memory
    """

    fastadict = defaultdict(dict)
    with open(ancAllele, 'r') as fa:
        for line in fa:
            x = line.split()
            fastadict[x[0]][int(x[1])] = x[2]

    fasta_sequences = SeqIO.parse(fastaFile, 'fasta')
    with open("ancRef.fasta", 'w') as out_file:
        for fasta in fasta_sequences:
            # read in header and sequence
            header, sequence = fasta.id, str(fasta.seq)
            seq = list(sequence)  # strings are immutable
            for pos in fastadict[header]:
                allele = fastadict[header][pos]
                # replace base w/ allele at pos-1
                if seq[pos-1].islower():
                    allele = allele.lower()
                seq[pos-1] = allele
            # when done with the header, write
            out_file.write(">{}\n{}\n".format(header, ''.join(seq)))
    return(None)


if __name__ == "__main__":
    vcf2fasta(args.fasta, args.ancestral)


