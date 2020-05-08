#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 09:58:25 2020
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
    This is an example of an indented section. It's like any other section,
    but the body is indented to help it stand out from surrounding text.

If a section is indented, then a section break is created by
resuming unindented text.

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
import re
import numpy as np
import glob
import argparse
import sys
import os


def buildMaskFile(infile, window, nLength, skipMask, numb):
    """Create mask from fasta or resample existing mask.

    Parameters
    ----------
    mask_name : TYPE
        DESCRIPTION.
    window : TYPE
        DESCRIPTION.
    nLength : TYPE
        DESCRIPTION.
    skipMask : TYPE
        DESCRIPTION.
    numb : TYPE
        DESCRIPTION.

    Returns
    -------
    mask_file : str
        mask file name.

    """
    mask_abspath = os.path.abspath(infile)
    mask_path, mask_name = os.path.split(mask_abspath)
    if "mask" in mask_name.split(".")[-1]:
        # need larger number for masking longer training sims
        mask_files = glob.glob(f"{os.path.join(mask_path, mask_name)}.*.mask")
        if mask_files:
            numb_list = []
            for file in mask_files:
                numb_list.append(int(file.split(".")[-2]))
            numb_max = max(numb_list)
            resample_file = f"{mask_name}.{numb_max}.mask"
            if numb_max < numb:
                with open(os.path.join(mask_path, resample_file), 'r') as mask:
                    coord_list = []
                    for line in mask:
                        if line.strip():
                            if not line.startswith("//"):
                                tmp_coord = []
                                while line.strip():
                                    tmp_coord.append(line)
                                    line = next(mask)
                                coord_list.append(tmp_coord)
                sample_file = f"{mask_name}.{numb}.mask"
                with open(os.path.join(mask_path, sample_file), 'w') as mask_renumb:
                    n = 0
                    for coord in coord_list:
                        for c in coord:
                            mask_renumb.write(c)
                        mask_renumb.write("\n//\n\n")
                        n += 1
                    for rand_coord in np.random.choice(coord_list, numb-n+1):
                        for coord in rand_coord:
                            for c in coord:
                                mask_renumb.write(c)
                        mask_renumb.write("\n//\n\n")
            else:
                return os.path.join(mask_path, resample_file)
    elif "fa" in mask_name.split(".")[-1]:
        fasta_file = f"{os.path.join(mask_path, mask_name)}"
        rN = re.compile(r'N{{{0},}}'.format(nLength))
        with open(f"{fasta_file}", 'r') as fasta:
            seq_list = []
            for line in fasta:
                if not line.startswith(">"):
                    seq_list.append(line.strip())
            seq = "".join(seq_list)
        n = 0
        sample_file = f"{mask_name}.{numb}.mask"
        with open(os.path.join(mask_path, sample_file), 'w') as mask_numb:
            # seqRlist = []
            coord_list = []
            # start mask
            start = 0
            end = start + window
            step = window
            while end < len(seq):
                seqR = seq[start:end]
                if (seqR.count('N') / window) <= skipMask:
                    coord = [(m.start(), m.end()) for m in re.finditer(rN, seqR)]
                    if coord:
                        coord_list.append(coord)
                        for i, j in coord:
                            mask_numb.write(f"0 {i/window} {j/window}\n")
                    else:
                        mask_numb.write("0 0 0\n")
                        coord_list.append([0, 0])
                    mask_numb.write("\n//\n\n")
                    # seqRlist.append(seqR)
                start += step
                end += step
                n += 1
            if n < numb:
                for rand_coord in np.random.choice(coord_list, numb-n+1):
                    for i, j in rand_coord:
                        mask_numb.write(f"0 {i/window} {j/window}\n")
                    mask_numb.write("\n//\n\n")
    return os.path.join(mask_path, sample_file)


def makeUnmaskedFrac(mask_file, file=True):
    """List of unmasked fraction for each windows.

    Parameters
    ----------
    mask_file : TYPE
        DESCRIPTION.

    Returns
    -------
    frac_mask : TYPE
        DESCRIPTION.

    """
    frac_mask = []
    with open(mask_file, 'r') as mask:
        for line in mask:
            if line.strip():
                if not line.startswith("//"):
                    frac = 0
                    while line.strip():
                        x = line.split()
                        frac += float(x[2]) - float(x[1])
                        line = next(mask)
                    frac_mask.append(1-frac)
    if file:
        path, mask_name = os.path.split(mask_file)
        with open(os.path.join(path, f"{mask_name}-unfracmask"), 'w') as out:
            for frac in frac_mask:
                out.write(f"{frac}\n")
    return frac_mask


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', "--infile", required=True, help="fasta/mask input file")
    parser.add_argument('-w', "--window_size", type=int, default=10000,
                        help="window size for masking fractions")
    parser.add_argument('-n', "--nLength", type=int, default=100,
                        help="length of Ns to cosider for masking")
    parser.add_argument('-m', "--skipMasking", type=float, default=0.50,
                        help="skip sites with more that this percent of sites as N")
    parser.add_argument('-nMask', "--numberMask", type=int,
                        help="how many mask lines to resample")
    return parser.parse_args(args_in)


def main():
    """Run main function."""
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    mask_file = args.infile
    window_size = args.window_size
    n_length = args.nLength
    skip_mask = args.skipMasking
    n_mask_lines = args.numberMask
    # =========================================================================
    #  Main executions
    # =========================================================================
    mask_filename = buildMaskFile(mask_file, window_size, n_length, skip_mask, n_mask_lines)
    frac_mask = makeUnmaskedFrac(mask_filename)


if __name__ == "__main__":
    main()
