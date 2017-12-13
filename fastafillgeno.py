#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 12:53:50 2017
@author: stsmall
vcf2geno #simont martin
fastaFillgeno
geno2vcf
"""
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-i', "--infile", type=str, required=True,
                    help="infile in genoformat")
parser.add_argument('-o', "--outfile", type=str, required=True,
                    help="outfile for filled geno")
parser.add_argument('-b', "--bed", type=str, required=True,
                    help="bed file from a fasta file, 1 base per site. Make"
                    "this by first making a bed file of all positions and then"
                    "using bedtools getFastafrombed -tab -Outbed")
parser.add_argument('-v', "--mafvcf", type=str, required=False,
                    help="optional vcf from maffilter for the outgroup fill"
                    "ensure that this is coord sorted since maffilter does not"
                    "do this by default")
parser.add_argument('-r', "--outgroup", nargs='+', required=False,
                    help="outgroup sample name(s)")
args = parser.parse_args()


def readbedfasta(fasta):
    """
    """
    fdict = defaultdict(dict)
    with open(fasta, 'r') as fbed:
        for line in fbed:
            x = line.strip().split()
            fdict[x[0]][x[2]] = x[-1]
    return(fdict)


def readmafvcf(vcf):
    """
    """
    vcfdict = defaultdict(dict)
    with open(vcf, 'r') as vbed:
        for line in vbed:
            if not line.startswith("#"):
                x = line.strip().split()
                vcfdict[x[0]][x[1]] = (x[3], x[4])
    return(vcfdict)


def fillgeno(infile, outfile, fdict, vcfdict, outgroup, dlm="."):
    """
    """
    fillposi = 1
    refcount = 0
    with open(outfile, 'w') as f:
        with open(infile, 'r') as geno:
            for line in geno:
                if line.startswith("#CHROM"):
                    f.write(line)
                    header = line.strip().split()
                    samples = header[2:]
                    n_samples = len(samples)
                    out_iix = []
                    for o in outgroup:
                        out_iix = [i for i, x in enumerate(header) if o == x.split(dlm)[0]]
                else:
                    x = line.strip().split()
                    chrom = x[0]
                    pos = int(x[1])
                    while fillposi < pos:
                        fillpos = str(fillposi)
                        try:
                            nuc = fdict[chrom][fillpos]
                            lnuc = "{}/{}".format(nuc, nuc)
                            reflist = [lnuc] * n_samples
                            if outgroup and vcfdict:
                                try:
                                    anc = vcfdict[chrom][fillpos][1]
                                    # ref = vcfdict[chrom][x[1]][0]
                                    for ix in out_iix:
                                        reflist[ix] = "{}/{}".format(anc, anc)
                                    refcount += 1
                                except KeyError:
                                    for ix in out_iix:
                                        import ipdb;ipdb.set_trace()
                                        if "/" in reflist[0]:
                                            reflist[ix] = "N/N"
                                        else:
                                            reflist[ix] = "N"
                            else:
                                pass
                            f.write("{}\t{}\t{}\n".format(chrom, fillpos,
                                    "\t".join(reflist)))
                        except KeyError:
                            continue
                        fillposi += 1
                    f.write(line)
    print("filled ref position in {}: {}".format(outgroup, refcount))
    return(None)


if __name__ == "__main__":
    infile = args.infile
    outfile = args.outfile
    fasta = args.bed
    vcf = args.mafvcf
    outgroup = args.outgroup
    fdict = readbedfasta(fasta)
    if args.mafvcf:
        vcfdict = readmafvcf(args.mafvcf)
    else:
        vcfdict = ''
    fillgeno(infile, outfile, fdict, vcfdict, outgroup)
