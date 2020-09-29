#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 14:17:01 2020
@author: Scott T. Small

This module demonstrates documentation as specified by the `NumPy
Documentation HOWTO`_. Docstrings may extend over multiple lines. Sections
are created with a section header followed by an underline of equal length.

Example
-------
Examples can be given using either the ``Example`` or ``Examples``
sections. Sections support any reStructuredText formatting, including
literal blocks::

    $ python example_numpy.py


Section breaks are created with two blank lines. Section breaks are also
implicitly created anytime a new section starts. Section bodies *may* be
indented:

Notes
-----
	count accessible w/ `*`preds coordinates in fasta
		>50% of the sites are N, then the site is skipped
		need to pull in 2 files, take the preds bed for the same file,
        move site by site only counting those with <50% Ns
		python slide_accessible.py 2L.preds 2L.van.fa 2L.fun.fa 0.50
			# preds :: get coords
			# read in FASTA with biopython
			# use 0.50 to check the % missing in the case it is not 50%

Attributes
----------
module_level_variable1 : int
    Module level variables may be documented in either the ``Attributes``
    section of the module docstring, or in an inline docstring immediately
    following the variable.

    Either form is acceptable, but the two should not be mixed. Choose
    one convention to document module level variables and be consistent
    with it.

"""
import sys
import argparse
from Bio import AlignIO
import numpy as np
from os import path
from tqdm import tqdm


def count_aln_N(count_N, records):
    """Test whether a specific base is an N or n. Return the count for a window.

    Parameters
    ----------
    count_N : TYPE
        DESCRIPTION.
    fasta1_aln : TYPE
        DESCRIPTION.
    fasta2_aln : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    for record in records:
        for i, seq in enumerate(record):
            if seq == "N" or seq == "n":
                count_N[i] += 1

    return count_N


def open_fasta(chrom, fasta1, fasta2):
    fasta1_dir, fasta1_base = path.split(fasta1)
    fasta2_dir, fasta2_base = path.split(fasta2)
    new_fasta1 = path.join(fasta1_dir, f"{chrom}.{fasta1_base}")
    new_fasta2 = path.join(fasta2_dir, f"{chrom}.{fasta2_base}")
    fasta1_aln = AlignIO.read(new_fasta1, 'fasta')
    len_f1 = len(fasta1_aln)
    fasta2_aln = AlignIO.read(new_fasta2, 'fasta')
    len_f2 = len(fasta2_aln)
    total_aln = len_f1 + len_f2
    assert fasta1_aln.get_alignment_length() == fasta2_aln.get_alignment_length()

    return fasta1_aln, fasta2_aln, total_aln


def count_accessible(fasta1, fasta2, start, end, total_aln, miss):
    window_size = end - start
    count_N = np.zeros(window_size)
    count_N = count_aln_N(count_N, fasta1[:, start:end])
    count_N = count_aln_N(count_N, fasta2[:, start:end])
    miss_sites = len(np.where(count_N/total_aln >= miss)[0])
    access_bp = window_size - miss_sites

    return access_bp


def find_accessible(preds, fasta1, fasta2, miss):
    """Calculate the number of accessible bases for FILET analysis.

    Parameters
    ----------
    preds : str
        File DESCRIPTION.
    fasta1 : str
        File DESCRIPTION.
    fasta2 : str
        File DESCRIPTION.
    miss : float
        DESCRIPTION.

    Returns
    -------
    None.

    """
    out = open("accessible.txt", 'w')
    count_line = 0
    with open(preds, 'r') as coords:
        line = coords.readline()
        chrom, start, end, med, prob = line.split()
        count_line += 1
        for line in coords:
            count_line += 1

    pbar = tqdm(total=count_line)
    with open(f"{preds}-accessible", 'w') as bed:
        f1, f2, total_aln = open_fasta(chrom, fasta1, fasta2)
        with open(preds, 'r') as coords:
            basepairs = 0
            for line in coords:
                pbar.update(1)
                new_chrom = line.split()[0]
                if chrom == new_chrom:
                    chrom, start, end, med, prob = line.split()
                    try:
                        access_bp = count_accessible(f1, f2, int(start), int(end), total_aln, miss)
                        basepairs += access_bp
                        bed.write(f"{chrom}\t{start}\t{end}\t{access_bp}\t{prob}\n")
                    except IndexError:
                        end = total_aln
                        access_bp = count_accessible(f1, f2, int(start), int(end), total_aln, miss)
                        basepairs += access_bp
                        bed.write(f"{chrom}\t{start}\t{end}\t{access_bp}\t{prob}\n")
                else:
                    out.write(f"{preds}\t{chrom}\t{basepairs}\n")
                    f1, f2, total_aln = open_fasta(new_chrom, fasta1, fasta2)
                    basepairs = 0
                    # catch that count
                    chrom, start, end, med, prob = line.split()
                    access_bp = count_accessible(f1, f2, int(start), int(end), total_aln, miss)
                    basepairs += access_bp
                    bed.write(f"{chrom}\t{start}\t{end}\t{access_bp}\t{prob}\n")

    out.write(f"{preds}\t{chrom}\t{basepairs}\n")
    pbar.close()
    out.close()
    return None


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--preds", type=str, help="preds or bed file with "
                        " coordinates")
    parser.add_argument("-f1", "--fasta1", type=str, help="fasta for pop1")
    parser.add_argument("-f2", "--fasta2", type=str, help="fasta for pop2")
    parser.add_argument("-m", "--siteMissFlt", type=float, help="percent missing "
                        "to not count site")
    return(parser.parse_args(args_in))


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    preds = args.preds
    fasta1 = args.fasta1
    fasta2 = args.fasta2
    miss = args.siteMissFlt
    # =========================================================================
    #  Main executions
    # =========================================================================
    find_accessible(preds, fasta1, fasta2, miss)


if __name__ == "__main__":
    main()
