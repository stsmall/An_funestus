#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 17:24:58 2019
parseTwisst.py -i dist.txt -t [topos] -p [pairs]

@author: stmall
"""
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', "--infile", type=str, required=True,
                    help="infile from distmatrix twisst")
parser.add_argument('-t', "--topos", nargs='+', required=False,
                    help="list of topos e.g., topo82")
parser.add_argument('-p', "--pairs", nargs='+', required=False,
                    help="list of interested pairs, e.g. Fun-Lik")
parser.add_argument('-f', "--minfreq", type=float, default=0.0)
parser.add_argument("--div", action="store_true", help="a divergence file")
parser.add_argument("--blen", action="store_true", help="a file of branch"
                    "lengths")
parser.add_argument("--tplot", action="store_true", help="boxplot of topos")
parser.add_argument("--pplot", action="store_true", help="boxplot of pairs")
parser.add_argument("--bplot", action="store_true", help="boxplot of nodeage")
args = parser.parse_args()


def boxplotD(div_df, pair, topo, topoplot, pairsplot):
    """
    """
    import ipdb;ipdb.set_trace()
    if pairsplot:
        label = pair
        x = [l for l in div_df.loc[pair]]
        poplist = list(div_df.loc[pair].index.values)
    elif topoplot:
        label = topo
        x = [l for l in div_df[topo]]
        poplist = list(div_df[topo].index.values)
    fig, ax = plt.subplots(figsize=(4, 4))
    sns.despine(ax=ax, offset=5)
    lw = 1.5
    # import ipdb; ipdb.set_trace()
    box = ax.boxplot(
        x,
        labels=poplist,  patch_artist=True,
        medianprops={"color": "k", "linewidth": lw},
        whiskerprops={"color": "k"},
        capprops={"color": "k"},
        showfliers=False,
        flierprops={"c": "k", "markersize": 2})
    ax.set_ylabel(r'Divergence', rotation=0, fontsize=16)
    # set list of colors
    colornames = list(mpl.colors.cnames.keys())[:len(poplist)]
    for patch, color in zip(box['boxes'], colornames):
        patch.set_facecolor(color)
        patch.set_linewidth(lw)
    fig.savefig("{}.pdf".format(label), bbox_inches="tight")
    return(None)


def boxplotB(data):
    """
    """
    poplist = ["topo{}".format(i) for i in range(1, 106)]
    fig, ax = plt.subplots(figsize=(4, 4))
    sns.despine(ax=ax, offset=5)
    lw = 1.5
    box = ax.boxplot(
        data,
        labels=poplist,  patch_artist=True,
        medianprops={"color": "k", "linewidth": lw},
        whiskerprops={"color": "k"},
        capprops={"color": "k"},
        showfliers=False,
        flierprops={"c": "k", "markersize": 2})
    ax.set_ylabel(r'Divergence', rotation=0, fontsize=16)
    # set list of colors
    colornames = list(mpl.colors.cnames.keys())[:len(poplist)]
    for patch, color in zip(box['boxes'], colornames):
        patch.set_facecolor(color)
        patch.set_linewidth(lw)
    fig.savefig("NodeDepth.pdf", bbox_inches="tight")
    return(None)


def getDivergence(infile, topos, pairs, toposplot, pairsplot):
    """parses out distance file to return a distribution for each species pair
    """
    div_dict = defaultdict(lambda: defaultdict(list))
    with open(infile,'r') as f:
        header = next(f).split()
        for line in f:
            for i, div in enumerate(line.split()):
                topo, ind1, ind2 = header[i].split("_")
                pairname = "{}-{}".format(ind1, ind2)
                div_dict[topo][pairname].append(float(div))
    div_df = pd.DataFrame(div_dict)
    
    # make boxplot
    if toposplot:
        for topo in topos:  # for each topo returns all pairwise distances on 1 plot
            boxplotD(div_df, pairs, topo, toposplot, pairsplot)
    if pairsplot:
        for pair in pairs:  # for each pair returns a boxplot of distances on each topo
            boxplotD(div_df, pair, topos, toposplot, pairsplot)
    
    # calculate mean distances
    f = open("meandist.out", 'w')
    if topos:  # specific topos
        for t in topos:
            if pairs:  # specific pairs
                for p in pairs:
                    pmean = np.nanmean(div_df[t].loc[p])
                    pNF = np.nanpercentile(div_df[t].loc[p], 97.5)
                    pTF = np.nanpercentile(div_df[t].loc[p], 2.5)
                    f.write("{}:{}:{} [{}-{}]\n".format(t, p, pmean, pTF, pNF))
            else:
                for p in list(div_df.index.values):
                    pmean = np.nanmean(div_df[t].loc[p])
                    pNF = np.nanpercentile(div_df[t].loc[p], 97.5)
                    pTF = np.nanpercentile(div_df[t].loc[p], 2.5)
                    f.write("{}:{}:{} [{}-{}]\n".format(t, p, pmean, pTF, pNF))
    else:
        for t in list(div_df.columns.values):
            if pairs:  # specific pairs
                for p in pairs:
                    pmean = np.nanmean(div_df[t].loc[p])
                    pNF = np.nanpercentile(div_df[t].loc[p], 97.5)
                    pTF = np.nanpercentile(div_df[t].loc[p], 2.5)
                    f.write("{}:{}:{} [{}-{}]\n".format(t, p, pmean, pTF, pNF))
            else:
                for p in list(div_df.index.values):
                    pmean = np.nanmean(div_df[t].loc[p])
                    pNF = np.nanpercentile(div_df[t].loc[p], 97.5)
                    pTF = np.nanpercentile(div_df[t].loc[p], 2.5)
                    f.write("{}:{}:{} [{}-{}]\n".format(t, p, pmean, pTF, pNF))
    f.close()
    return(None)


def sumBranchLengths(infile, minFreq, nodedepthplot, step=10, topos=105):
    """
    """
    blen_box = []
    tree_count = 0
    with open(infile, 'r') as f:
        for line in f:
            if line.startswith("#"):
                pass
            else:
                tree_count += 1
                i = 0
                j = step
                blen = list(map(float, line.split()))
                blen_list = []
                while j <= len(blen):
                    blen_list.append(sum(blen[i+1:j]))  # this skips the first which is outgroup BL
                    i += step
                    j += step
                blen_box.append(blen_list)
    
    # boxplot of node depth
    data = list(zip(*blen_box))
    
    if nodedepthplot:
        boxplotB(data)
    
    # output file
    f = open("nodedepth.out", 'w')
    for t in range(1, 106):
        nancount = sum(np.isnan(data[t-1]))
        if (1 - (nancount/tree_count)) >= minFreq:
            m = np.nanmean(data[t-1])
            mn = np.nanmedian(data[t-1])
            pl = np.nanpercentile(data[t-1], 2.5)
            pu = np.nanpercentile(data[t-1], 97.5)
            f.write("topo{}:{} {} [{}-{}]\n".format(t, m, mn, pl, pu))
    f.close()  
    
    return(None)


if __name__ == "__main__":
    infile = args.infile
    topos = args.topos
    pairs = args.pairs
    if args.div:
        getDivergence(infile, topos, pairs, toposplot=args.tplot, pairsplot=args.pplot)
    elif args.blen:
        sumBranchLengths(infile, args.minfreq, nodedepthplot=args.bplot)
