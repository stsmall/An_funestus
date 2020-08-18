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
		need to pull in 2 files, take the preds bed for the same file, move site by site only counting those with <50% Ns
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
from tqdm.contrib import tenumerate


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
    # build coords list
    coord_list = []
    with open(preds, 'r') as coords:
        for line in coords:
            line = line.split()
            chrom = line[0]
            start = int(line[1])
            coord_list.append(start)
    # load alignments
    fasta1_aln = AlignIO.read(fasta1, 'fasta')
    len_f1 = len(fasta1_aln)
    fasta2_aln = AlignIO.read(fasta2, 'fasta')
    len_f2 = len(fasta2_aln)
    total_aln = len_f1 + len_f2
    assert fasta1_aln.get_alignment_length() == fasta2_aln.get_alignment_length()

    # count N's
    count_N = np.zeros([total_aln, fasta1_aln.get_alignment_length()], dtype=object)
    i = 0
    for record in fasta1_aln:
        count_N[i] = np.array(list(record.seq), dtype=object) == "N"
        i += 1
    for record in fasta2_aln:
        count_N[i] = np.array(list(record.seq), dtype=object) == "N"
        i += 1
    breakpoint()
    # mask_N = count_N == "N"
    sum_N = np.sum(count_N, axis=1)

    # count N's per window
    basepairs = 0
    with open(f"{preds}-accessible", 'w') as bed:
        for i, start in tenumerate(coord_list):
            try:
                end = coord_list[i+1]
                window_size = end - start
                miss_count = sum_N[start:end]
                miss_sites = len(np.where(miss_count/total_aln >= miss)[0])
                access_bp = window_size - miss_sites
                basepairs += access_bp
                bed.write(f"{chrom}\t{start}\t{end}\t{access_bp}\n")
            except IndexError:
                break
    print(f"{chrom}: {basepairs}")


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--preds", type=str, help="preds or bed file with "
                        " coordinates")
    parser.add_argument("-f1", "--fasta1", type=str, help="fasta for pop1")
    parser.add_argument("-f2", "--fasta2", type=str, help="fasta for pop2")
    parser.add_argument("-m", "--siteMissFlt", type=float, default=0.50,
                        help="percent missing to not count site")
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
