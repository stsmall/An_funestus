#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 11:42:48 2018
fillMAF.py MAF vcf.list

suggestions:
    hal2maf
    filterDuplicates
    maffilter: subset, concant, etc ...
    fillMAF.py
    maffilter: windows
@author: stsmall
"""

import argparse
from collections import defaultdict
import re
import random
import ipdb

parser = argparse.ArgumentParser()

parser.add_argument('-v', "--vcfList", type=str, required=True,
                    help='List of VCF files')
parser.add_argument('-m', "--mafFile", required=True, type=str,
                    help='path to maf file')
parser.add_argument('-U', '--iupac', action='store_true',
                    help='add IUPAC base for het')
parser.add_argument('-M', '--major', action='store_true',
                    help='call major allel if het based on coverage')
parser.add_argument('-R', '--rand', action='store_true',
                    help='choose random allele at hets')
parser.add_argument('-N', '--Nasmissing', action='store_true',
                    help='add N as missing, default is ref')
args = parser.parse_args()


def transCoord(coord):
    """Translates maf coordinates to VCF format
    """
    p = re.compile("N")
    start, aln, strand, length, seq = coord
    length = int(length)
    pos_list = []
    poslist = []
    # find all "N", return coords as list
    seq1 = seq.replace("-", "")
    for m in p.finditer(seq1):
        pos_list.append(m.start())  # in vcf coords
    for m in p.finditer(seq):
        poslist.append(m.start())  # in maf coords
    if "-" in strand:
        n_list = [length - i for i in pos_list]
    else:
        n_list = [i + 1 for i in pos_list]
    return(n_list, poslist)


def getMAFambig(mafFile):
    """Find all N positions in mafFile

    Parameters
    ------
    mafFile: file, maf formated file

    Returns
    ------
    mafdict: dict, dictionary storing N positions in VCF, 1based corrdinates
        [IND][CHROM]: (Key, coord, "N")

    """
    mafdict = defaultdict(lambda: defaultdict(lambda: []))
    with open(mafFile, 'r') as maf:
        for line in maf:
            try:
                if line.startswith("a"):
                    line = next(maf)
                    while line.startswith("s"):
                        x = line.split()
                        ind, chrom = x[1].split(".")[:2]
                        block_key = "{}".format("_".join(x[1:6]))
                        coord, poslist = transCoord(x[2:])  # coord is list of N in block
                        if coord:
                            mafdict[ind][chrom].append([block_key, coord, poslist, "N"])
                        line = next(maf)
            except StopIteration:
                break
    return(mafdict)


def IUPAC(ALLELE):
    """ALLELE needs to be transformed into IUPAC
    """
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
        ALLELE = "X"
    return(ALLELE)


def getVCFnucleotide(vcfFile, iupac, rand, major, Nasmissing):
    """Get allele from VCF
    """
    vcfdict = defaultdict(dict)
    with open(vcfFile, 'r') as vcfList:
        for vcf in vcfList:
            with open(vcf, 'r') as v:
                for line in v:
                    if not line.startswith("##"):
                        if line.startswith("#CHROM"):
                            samples = line.split()
                        else:
                            x = line.split()
                            for i, s in enumerate(x[9:]):
                                gt = s.split(":")
                                if gt[0].count("1") == 2:
                                    allele = x[4]
                                elif gt[0].count("0") == 2:
                                    allele = x[3]
                                else:
                                    allele = x[3] + x[4]
                                if len(allele) > 1:
                                    if iupac:
                                        allele = IUPAC(allele)
                                    elif rand:
                                        allele = random.choice(allele)
                                    elif major:
                                        try:
                                            # find AD field
                                            af = map(int, gt[1].split(","))
                                            # max of AD
                                            af_ix = af.index(max(af))
                                            # stupid tri-allelic
                                            if len(af) > 2:
                                                if af_ix == 0:
                                                    allele = x[3]
                                                elif af_ix == 1:
                                                    allele = x[4].split(",")[0]
                                                else:
                                                    allele = x[4].split(",")[1]
                                            else:
                                                allele = x[af_ix + 3]
                                        except ValueError:
                                            if Nasmissing:
                                                allele = "N"
                                            else:
                                                allele = x[3]
                                indname = samples[i+9].replace(".", "_")
                                vcfdict[indname]["{}_{}".format(x[0], x[1])] = allele
    return(vcfdict)


def replaceMaf(vcfdict, mafdict):
    """Replace "N" in mafdict with allele from vcfdict
    """
    for ind in mafdict.keys():
        for chrom in mafdict[ind].keys():
            for i, coord in enumerate(mafdict[ind][chrom]):
                pos = coord[1]
                a = ''
                for p in pos:
                    a += vcfdict[ind][chrom+"_"+pos]
                mafdict[ind][chrom][0][i][-1] = a
    return(mafdict)


def fillMaf(mafdict, mafFile):
    """Fill Ns in maf file

    mafdict[ind][chrom].append([block_key, coord, "N"])

    """
    f = open("{}.fill".format(mafFile), 'w')
    with open(mafFile, 'r') as maf:
        for line in maf:
            if line.startswith("a"):
                f.write("a\n")
                line = next(maf)
                while line.startswith("s"):
                    x = line.split()
                    ind, chrom = x[1].split(".")
                    seq = list(x[-1])
                    # find block key in mafdict list
                    block_key = x[1]+"_"+x[2]+"_"+x[3]+"_"+x[4]+"_"+x[5]
                    blocklist = zip(*mafdict[ind][chrom][0])[0]
                    i = blocklist.index(block_key)
                    coord = mafdict[ind][chrom][0][i]
                    for pos, nuc in zip(coord[2], coord[-1]):
                        seq[pos] = nuc
                    x[-1] = "".join(seq)
                    f.write("{}\n".format("\t".join(x)))
                    line = next(maf)
                f.write("\n")
    f.close()
    return(None)


if __name__ == "__main__":
    mafdict = getMAFambig(args.mafFile)
    vcfdict = getVCFnucleotide(args.vcfList, args.iupac, args.rand, args.major, args.Nasmissing)
    mafdict = replaceMaf(vcfdict, mafdict)
    fillMaf(mafdict, args.mafFile)
