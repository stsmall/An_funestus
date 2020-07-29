#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:24:44 2020
@author: Scott T. Small

This program creates a consensus fasta sequence, new reference, that uses inform-
ation from a population vcf to replace the reference base with fixed or high freq
ALT calls

Example
-------
Examples can be given using either the ``Example`` or ``Examples``
sections. Sections support any reStructuredText formatting, including
literal blocks::

    $ python example_numpy.py


Section breaks are created with two blank lines. Section breaks are also
implicitly created anytime a new section starts. Section bodies *may* be
indented:


"""
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict
import gzip
import re


def IUPAC(ALLELE):
    """ALLELE needs to be transformed into IUPAC."""
    # IUPAC transformation
    if "A" in ALLELE and "G" in ALLELE:
        ALLELE = "R"
    elif "C" in ALLELE and "T" in ALLELE:
        ALLELE = "Y"
    elif "G" in ALLELE and "C" in ALLELE:
        ALLELE = "S"
    elif "A" in ALLELE and "T" in ALLELE:
        ALLELE = "W"
    elif "G" in ALLELE and "T" in ALLELE:
        ALLELE = "K"
    elif "A" in ALLELE and "C" in ALLELE:
        ALLELE = "M"
    else:
        ALLELE = "N"

    return ALLELE


def dict_to_fasta(vcfdict, fasta_ref):
    """Read in fasta, changes base using dictionary entry.

    note that fasta will be 0 based so position wll be -1
    notes on Biopython: SeqIO.to_dict() which builds all sequences into a
    dictionary and save it in memory
    SeqIO.index() which builds a dictionary without putting the sequences
    in memory

    Parameters
    ----------
    vcfdict : TYPE
        DESCRIPTION.
    fasta_ref : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    fasta_sequences = SeqIO.parse(fasta_ref, 'fasta')
    with open(f'{fasta_ref}-consensus', 'w') as out_file:
        for fasta in fasta_sequences:
            chrom, sequence = fasta.id, str(fasta.seq)
            seq = list(sequence)
            for pos in vcfdict[chrom].keys():
                new_allele = vcfdict[chrom][pos][1]
                ref = vcfdict[chrom][pos][0]
                if seq[pos-1] != "N":
                    if seq[pos-1].islower():
                        assert seq[pos-1] == ref.lower()
                        new_allele = new_allele.lower()
                    assert seq[pos-1] == ref
                    seq[pos-1] = new_allele
            out_file.write(f">{chrom}\n{''.join(seq)}\n")


def freq_to_dict(frq_file, aaf, iupac, ploidy=2):
    """Read a --freq file from vcftools to a dictionary.

    Parameters
    ----------
    frq_file : TYPE
        DESCRIPTION.
    aaf : TYPE
        DESCRIPTION.
    iupac : TYPE
        DESCRIPTION.
    ploidy : TYPE, optional
        DESCRIPTION. The default is 2.

    Returns
    -------
    vcfdict : TYPE
        DESCRIPTION.

    """
    frqdict = defaultdict(dict)

    if frq_file.endswith(".gz"):
        fopen = gzip.open
    else:
        fopen = open

    with fopen(frq_file, 'rt') as frq:
        for line in frq:
            if line.startswith("CHROM"):
                header = line.strip().split()
            else:
                frqline = line.strip().split()
                chrom = frqline[0]
                pos = int(frqline[1])
                nalleles = frqline[2]
                nchroms = frqline[3]
                ref, fref = frqline[4].split(":")
                alt, falt = frqline[5].split(":")
                if "nan" in fref or "nan" in falt:
                    continue
                else:
                    falt = float(falt)
                    if len(frqline) == 7:
                        alt2, falt2 = frqline[6].split(":")
                        falt2 = float(falt2)
                        if falt2 > falt:
                            alt = alt2
                            falt = falt2
                    elif len(frqline) == 8:
                        alt2, falt2 = frqline[6].split(":")
                        falt2 = float(falt2)
                        alt3, falt3 = frqline[7].split(":")
                        falt3 = float(falt3)
                        if falt3 > falt2 > falt:
                            alt = alt3
                            falt = falt3
                        elif falt2 > falt3 > falt:
                            alt = alt2
                            falt = falt2
                    if float(falt) >= aaf:
                        new_allele = alt
                        frqdict[chrom][pos] = ref + new_allele
                    elif iupac:
                        new_allele = IUPAC(ref + alt)
                        frqdict[chrom][pos] = ref + new_allele

    return frqdict


def vcf_to_dict(vcf_file, aaf, iupac, ploidy=2):
    """Read a vcf file and stores info in a dictionary.

    Here using a tuple for the key

    Parameters
    ----------
    vcf_file : TYPE
        DESCRIPTION.
    aaf : TYPE
        DESCRIPTION.
    iupac : TYPE
        DESCRIPTION.
    ploidy : TYPE, optional
        DESCRIPTION. The default is 2.

    Returns
    -------
    vcfdict : TYPE
        DESCRIPTION.
    """
    # TODO: add option to take list of samples
    vcfdict = defaultdict(dict)

    if vcf_file.endswith(".gz"):
        fopen = gzip.open
    else:
        fopen = open

    with fopen(vcf_file, 'rt') as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                samples = line.strip().split()
            if not line.startswith("#"):
                vcfline = line.strip().split()
                chrom = vcfline[0]
                pos = int(vcfline[1])
                ref = vcfline[3]
                alt = vcfline[4]
                haps = len(vcfline[9:]) * ploidy
                if "," in alt:
                    alt = vcfline[4].split(",")
                    freq1 = sum([re.split("/|\|", ind.split(":")[0]).count('1') for ind in vcfline[9:]])
                    freq2 = sum([re.split("/|\|", ind.split(":")[0]).count('2') for ind in vcfline[9:]])
                    if freq1 >= freq2:
                        alt = alt[0]
                        freq = freq1
                    else:
                        alt = alt[1]
                        freq = freq2
                else:
                    freq = sum([re.split("/|\|", ind.split(":")[0]).count('1') for ind in vcfline[9:]])

                if freq/haps >= aaf:
                    new_allele = alt
                    vcfdict[chrom][pos] = ref + new_allele
                elif iupac:
                    new_allele = IUPAC(ref + alt)
                    vcfdict[chrom][pos] = ref + new_allele

    return vcfdict


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcfFile', type=str, help='path to vcf file')
    parser.add_argument('--freqFile', type=str, help='path to freq file')
    parser.add_argument('--refFasta', type=str, help='path to fasta file')
    parser.add_argument('--alt_allele_freq', type=float, default=1.0,
                        help='ALT allele freq above which to consider fixed')
    parser.add_argument('--iupac', action='store_true',
                        help='add IUPAC base for het')
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    vcf_file = args.vcfFile
    freq_file = args.freqFile
    fasta_ref = args.refFasta
    aaf = args.alt_allele_freq
    iupac = args.iupac
    # =========================================================================
    #  Main executions
    # =========================================================================
    if vcf_file:
        seqdict = vcf_to_dict(vcf_file, aaf, iupac)
    elif freq_file:
        seqdict = freq_to_dict(freq_file, aaf, iupac)
    else:
        raise Exception("vcf or freq file must be provided!")

    dict_to_fasta(seqdict, fasta_ref)


if __name__ == "__main__":
    main()
