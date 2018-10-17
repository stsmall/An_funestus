#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 16:42:26 2017
show-coords -Trlcd
sort -k1,1 -k2,2n
@author: scott
"""
import argparse
from collections import defaultdict
# from collections import OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument("--infish", type=str, required=True, help="")
parser.add_argument("--outfile", type=str, required=True, help="to-from")
parser.add_argument("--infile", type=str, required=True, help="")
parser.add_argument("--region", action="store_true", help="")
parser.add_argument("--pacbio", action="store_true", help="")
parser.add_argument("--mapping", action="store_true", help="")
args = parser.parse_args()


def nonblank_lines(f):
    """
    """
    for l in f:
        line = l.rstrip()
        if line:
            yield line


def makeheader(h):
    """
    """
    header = [i.lstrip('[').rstrip(']') for i in h]
    return(header)


def makemap(infish, region, pacbio):
    """
    """
    ddfish = defaultdict(list)
    if not pacbio:
        with open(infish, 'r') as fish:
            for line in fish:
                x = line.strip().split(",")
                chrom = x[3].split(":")[0]
                sector = x[3].split(":")[1].split("-")[0]
                if region:
                    ddfish["{}_{}".format(chrom, sector)].append(x[1])
                else:
                    ddfish[chrom].append(x[1])
    else:
        with open(infish, 'r') as fish:
            for line in fish:
                x = line.strip().split()
                chrom = x[1].split("_")[0]
                sector = x[1].split("_")[1]
                if region:
                    ddfish["{}_{}".format(chrom, sector)].append(x[0])
                else:
                    ddfish[chrom].append(x[0])
    return(ddfish)


def addmaptoaln(infile, ddfish):
    """
    """
    with open("{}.map".format(infile), 'w') as f:
        with open(infile, 'r') as nuc:
            for line in nonblank_lines(nuc):
                if line.startswith('['):
                    f.write("{}\n".format(line))
                elif line.strip().split()[0].isdigit():
                    x = line.strip().split()
                    if any([x[-1] in k for k in ddfish.values()]):
                        for key in ddfish.keys():
                            if x[-1] in ddfish[key]:
                                chrom = "{}_{}".format(key, x[-1])
                                x[-1] = chrom
                                f.write("{}\n".format("\t".join(x)))
                                break
                    else:
                        f.write("{}\n".format(line))
                else:
                    f.write("{}\n".format(line))
    return(None)


def makelinks(ddfish, outfile, infile, size=1000):
    """links are from - to, but files are named to - from
    """
    qs = []
    qn = []
    ss = []
    sn = []
    with open("circos.{}.links.txt".format(outfile), 'w') as f:
        with open(infile, 'r') as nuc:
            for line in nonblank_lines(nuc):
                if line.startswith('['):
                    h = line.strip().split("\t")
                    print(h)
                    header = makeheader(h)
                elif line.strip().split()[0].isdigit():
                    x = line.strip().split()
                    alnlen = header.index("LEN 1")
                    if int(x[alnlen]) >= size:
                        if any([x[-1] in k for k in ddfish.values()]):
                            for key in ddfish.keys():
                                if x[-1] in ddfish[key]:
                                    chrom = "{}_{}".format(key, x[-1])
                                    x[-1] = chrom
                                    break
                        qn.append(x[-1])
                        lenq = header.index("LEN Q")
                        qs.append(x[lenq])
                        lens = header.index("LEN R")
                        ss.append(x[lens])
                        sn.append(x[-2])
                        rstart = header.index("S1")
                        rend = header.index("E1")
                        qstart = header.index("S2")
                        qend = header.index("E2")
                        import ipdb;ipdb.set_trace()
                        f.write("{}\n".format(
                                " ".join([x[-1], x[qstart],
                                          x[qend], x[-2],
                                          x[rstart],
                                          x[rend]])))
    return(zip(sorted(set(sn), key=sn.index), sorted(set(ss), key=ss.index)),
           zip(sorted(set(qn), key=qn.index), sorted(set(qs), key=qs.index)))


def makechrom(outfile, aln):
    """
    """
    for i, j in enumerate(outfile.split("-")):
        chromlist = aln[i]
        with open("circos.{}.{}.karyotype.txt".format(j, outfile), 'w') as f:
            for n, s in chromlist:
                f.write("chr - {} {} 0 {} black\n".format(n, n, s))


if __name__ == "__main__":
    fishin = args.infish
    outfile = args.outfile
    infile = args.infile
    ddfish = makemap(fishin, args.region, args.pacbio)
    addmaptoaln(infile, ddfish)
    if not args.mapping:
        aln = makelinks(ddfish, outfile, infile)
        makechrom(outfile, aln)
