#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 16:01:12 2018

@author: scott
"""
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile',
                    help='path to infile', required=True)
parser.add_argument('-p', "--phylonet", required=True,
                    help="phylonet command line")
parser.add_argument('--seq', action='store_true', help='MSMC_SEQ')
parser.add_argument("--gt", action='store_true', help="MCMC_GT")
args = parser.parse_args()


def parseTreelist(treeFile):
    """Parse newick with many trees to phylonet
    """
    treelist = []
    with open(treeFile, 'r') as tree:
        for line in tree:
            treelist.append(line.strip())
    return(treelist)


def phylonetTree(treelist, phylonetcmd):
    """
    """
    f = open("phylonet.GT.nex", 'w')
    f.write("#NEXUS\n\nBEGIN TREES;\n")
    for i, t in enumerate(treelist):
        f.write("Tree gt{} = {}\n".format(i, t))
    f.write("\nEND;\n\nBEGIN PHYLONET;\n")
    f.write("{}\n\nEND;".format(phylonetcmd))
    f.close()
    return(None)


def parseSeqFile(seqFile):
    """Parse Nexus formatted file
    """
    total_char = 0
    seqdict = defaultdict(list)
    with open(seqFile, 'r') as seq:
        for line in seq:
            if "dimensions" in line:
                _, tax, char = line.split()
                ntax = tax.split("=")[1]
                nchar = char.split("=")[1]
                total_char += int(nchar[:-1])
            elif line.startswith("matrix"):
                line = seq.next()
                locus = "{}_{}".format(line.split()[0].split(":")[1], nchar)
                try:
                    while not line.startswith("#NEXUS"):
                        seqdict[locus].append(line)
                        line = seq.next()
                except StopIteration:
                    break
    return(ntax, total_char, seqdict)


def phylonetSeq(ntax, tchar, seqdict, phylonetcmd):
    """
    """
    f = open("phylonet.SEQ.nex", 'w')
    f.write("#NEXUS\n")
    f.write("Begin data;\n\tDimensions ntax={} nchar={};\n".format(ntax, tchar))
    f.write('\tFormat datatype=dna symbols="ACTG" missing=? gap=-;\n')
    f.write("\tMATRIX\n")
    for locus in seqdict.keys():
        l, s = locus.split("_")
        f.write("[{}, {}]\n".format(l, s))
        for seq in seqdict[locus]:
            try:
                name, dna = seq.split()
            except ValueError:
                import ipdb;ipdb.set_trace()
            f.write("{} {}\n".format(name.split(":")[0], dna))
    f.write("\n;END;\n\nBEGIN PHYLONET;\n")
    f.write("{}\n\nEND;".format(phylonetcmd))
    f.close()
    return(None)


if __name__ == "__main__":
    if args.gt:
        treelist = parseTreelist(args.infile)
        phylonetTree(treelist, args.phylonet)
    elif args.seq:
        ntax, tchar, seqdict = parseSeqFile(args.infile)
        phylonetSeq(ntax, tchar, seqdict, args.phylonet)
