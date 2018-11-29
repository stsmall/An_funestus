#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 13:05:57 2018
BppScrape.py -p nCDS -s bpp.out -c 100000
returns: Coordinate File, Twisst-style weights file
@author: scott
"""
from __future__ import print_function
from __future__ import division
import argparse
import glob
from collections import defaultdict
import logging
import re

#logging.debug('This message should appear on the console')
#logging.info('So should this')
#logging.warning('And this, too')

parser = argparse.ArgumentParser()
parser.add_argument('-p', "--prefix", required=True, help="prefix of input"
                    "file")
parser.add_argument('-s', "--suffix", required=True, help="suffix of input"
                    "file")
parser.add_argument("--scafs", required=True, nargs='+', action="append",
                    help="scaffold or chromosome")
parser.add_argument('-c', "--chainLen", default=0, type=int, help="chain"
                    "length, default"
                    "or value of 0 will return proportions")
parser.add_argument("--loglevel", help="info, debug, warning, error, critical")
parser.add_argument("--log", action="store_false")
args = parser.parse_args()


def scrapeBpp(prefix, suffix, chain, scafs):
    """
    """
    import ipdb;ipdb.set_trace()
    topoList = []
    weightsDict = defaultdict(dict)
    for s in scafs:
        fileList = glob.glob("{}{}*{}".format(prefix, s, suffix))
        for bppOut in fileList:
            coord = bppOut.replace(prefix, '').replace(suffix, '')
            keycoord = s + ":" + coord
            with open(bppOut, 'r') as bpp:
                for line in bpp:
                    if line.startswith("(A)"):
                        for line in bpp:
                            if line.strip() == "":
                                break
                            else:
                                x = line.split()
                                topo = "".join(x[3:])
                                if chain > 0:
                                    weightsDict[keycoord][topo] = float(x[1]) * chain
                                else:
                                    weightsDict[keycoord][topo] = float(x[1])
                                topoList.append(topo)
                    else:
                        pass
    return(weightsDict, topoList)


def writeWeights(weightsDict, topoList):
    """
    """
    # topoX\t NEWICK\n
    # topo1\ttopo2\t ... \n
    # weights\t ... \n
    t = open("topos.out", 'w')
    f = open("weights.out", 'w')
    topoSet = list(set(topoList))
    topoHeader = ''
    topoCount = len(topoSet)
    for i, topo in enumerate(topoSet):
        t.write("#topo{}\t{}\n".format(i+1, topo))
        topoHeader += "topo{}\t".format(i+1)
    t.close()
    f.write(topoHeader + "\n")
    for coord in weightsDict.keys():
        scaf, start, stop = re.findall(r'\w+', coord)
        topolist = [0] * topoCount
        for topo in weightsDict[coord]:
            ix = topoSet.index(topo)
            topolist[ix] = weightsDict[coord][topo]
        f.write("{}\t{}\t{}\t{}\n".format(scaf, start, stop, "\t".join(map(str, topolist))))
    f.close()
    return(None)


if __name__ == "__main__":
    weightsDict, topoList = scrapeBpp(args.prefix, args.suffix, args.chainLen, args.scafs)
    writeWeights(args.scaf, weightsDict, topoList)
