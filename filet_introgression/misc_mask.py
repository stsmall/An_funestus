#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 22:12:46 2020
@author: Scott T. Small

"""
import re
import numpy as np
import glob


# if mask:
#     mask_name = f"{pop1}-{pop2}"
#     mask_filename = buildMaskFile(mask_name, 10000, 100, 0.50, 100)
#     frac_mask = makeUnmaskedFrac(mask_filename)




def buildMaskFile(mask_name, window, nLength, skipMask, numb):
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
    mask_files = glob.glob(f"{mask_name}.*.mask")
    if mask_files:
        numb_list = []
        for file in mask_files:
            numb_list.append(int(file.split(".")[-2]))
        numb_max = max(numb_list)
        resample_file = glob.glob(f"{mask_name}.{numb_max}.mask")
        if numb_max < numb:
            with open(resample_file, 'r') as mask:
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
            with open(sample_file, 'w') as mask_renumb:
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
            return sample_file
        else:
            return resample_file
    else:
        # make mask files
        rN = re.compile(r'N{{{0},}}'.format(nLength))
        with open(f"{mask_name}.fasta", 'r') as fasta:
            seq_list = []
            for line in fasta:
                if not line.startswith(">"):
                    seq_list.append(line.strip())
            seq = "".join(seq_list)
        n = 0
        sample_file = f"{mask_name}.{numb}.mask"
        with open(sample_file, 'w') as mask_numb:
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
    return sample_file


def makeUnmaskedFrac(mask_file):
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
    with open(mask_file, 'r') as m:
        for line in m:
            if line.strip():
                if not line.startswith("//"):
                    frac = 0
                    while line.strip():
                        x = line.split()
                        frac += float(x[2]) - float(x[1])
                        line = next(m)
                    frac_mask.append(1-frac)
    return frac_mask
