#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 16:42:26 2017

@author: scott
"""
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--infish", type=str, required=True, help="")
parser.add_argument("--outfile", type=str, required=True, help="to-from")
parser.add_argument("--infile", type=str, required=True, help="")
parser.add_argument("--region", action="store_true", help="")
args = parser.parse_args()


def makemap(infish, region):
    """
    """
    ddfish = defaultdict(list)
    try:
        with open(infish, 'r') as fish:
            for line in fish:
                x = line.strip().split(",")
                chrom = x[3].split(":")[0]
                sector = x[3].split(":")[1].split("-")[0]
                if region:
                    ddfish["{}_{}".format(chrom, sector)].append(x[1])
                else:
                    ddfish[chrom].append(x[1])
    except:
        import ipdb;ipdb.set_trace()
    return(ddfish)


def makelinks(ddfish, outfile, infile, size=5000):
    """links are from - to, but files are named to - from
    """
    qs = []
    qn = []
    ss = []
    sn = []
    with open("circos.{}.links.txt".format(outfile), 'w') as f:
        with open(infile, 'r') as nuc:
            for line in nuc:
                if line.strip().split()[0].isdigit():
                    x = line.strip().split()
                    if int(x[4]) >= size:
                        for key in ddfish.keys():
                            if x[11] in ddfish[key]:
                                chrom = "{}_{}".format(key, x[11])
                                x[11] = chrom
                                qn.append(chrom)
                                qs.append(x[7])
                                ss.append(x[6])
                                sn.append(x[10])
                                f.write("{}\n".format(" ".join(x[11], x[2],
                                                               x[3], x[10],
                                                               x[0], x[1])))
    return(zip(set(qn), set(qs)), zip(set(sn), set(ss)))


def makechrom(outfile, aln):
    """
    """
    for i, j in enumerate(outfile.split("-")):
        chromlist = aln[j]
        with open("circos.{}.{}.karyotype.txt".format(i, outfile), 'w') as f:
            for n, s in chromlist:
                f.write("chr - {} {} 0 {} black\n".format(n, n, s))


if __name__ == "__main__":
    fishin = args.infish
    outfile = args.outfile
    infile = args.infile
    aln = makelinks(makemap(fishin, args.region), outfile, infile)
    makechrom(outfile, aln)
