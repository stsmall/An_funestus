#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 17:24:58 2019
parseTwisst.py -i dist.txt -t [topos] -p [pairs]
will parse a divergence and branch length file produced by twisst
@author: stmall
"""
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
mpl.rcParams['pdf.fonttype'] = 42

parser = argparse.ArgumentParser()
parser.add_argument('-i', "--infile", type=str, required=True,
                    help="infile from distmatrix twisst")
parser.add_argument('-t', "--topos", nargs='+', required=False,
                    help="list of topos e.g., topo82")
parser.add_argument('-p', "--pairs", nargs='+', required=False,
                    help="list of interested pairs, e.g. Fun-Lik")
parser.add_argument('-f', "--minFreq", type=float, default=0.0)
parser.add_argument("--dist", action="store_true", help="a divergence file")
parser.add_argument("--blen", action="store_true", help="a file of branch"
                    "lengths")
parser.add_argument("--tplot", action="store_true", help="boxplot of topos")
parser.add_argument("--pplot", action="store_true", help="boxplot of pairs")
parser.add_argument("--bplot", action="store_true", help="boxplot of nodeage")
args = parser.parse_args()


def boxplotD(div_df, pair, topo, topoplot, pairsplot):
    """
    """
    if pairsplot:
        label = pair
        if topo:
            df_list = div_df[topo].loc[pair].sort_index()
            x = [l for l in df_list]
            poplist = list(df_list.index.values)
        else:
            df_list = div_df.loc[pair].sort_index()
            x = [l for l in df_list]
            poplist = list(df_list.index.values)            
    elif topoplot:
        label = topo
        x = [l for l in div_df[topo]]
        poplist = list(div_df[topo].index.values)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.despine(ax=ax, offset=5)
    lw = 1.5
    box = ax.boxplot(
        x,
        labels=poplist,  patch_artist=True,
        medianprops={"color": "k", "linewidth": lw},
        whiskerprops={"color": "k"},
        capprops={"color": "k"},
        showfliers=False,
        flierprops={"c": "k", "markersize": 2})
    ax.set_ylabel(r'PwDiv', rotation=90, fontsize=12)
    ax.set_xticklabels(poplist, rotation=45, fontsize=8)
    # set list of colors
    colornames = list(mpl.colors.cnames.keys())[:len(poplist)]
    for patch, color in zip(box['boxes'], colornames):
        patch.set_facecolor(color)
        patch.set_linewidth(lw)
    fig.savefig("{}.pdf".format(label), bbox_inches="tight")
    return(None)


def boxplotB(data, topos):
    """
    """
    if topos:
        poplist = topos
        topo_ix = [(int(x.strip("topo"))-1) for x in poplist]
        data_arr = np.array(data)[topo_ix]
        data_list = [a1[~np.isnan(a1)] for a1 in data_arr]
    else:
        poplist = ["topo{}".format(i) for i in range(1, 106)]
        data_arr = np.array(data)
        data_list = [a1[~np.isnan(a1)] for a1 in data_arr]
        good_ix = [i for i,j in enumerate(data_list) if len(j) > 0]
        poplist = np.array(poplist)[good_ix]
        data_list = [j for j in data_list if len(j) > 0]
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.despine(ax=ax, offset=5)
    lw = 1.5
    box = ax.boxplot(
        data_list,
        labels=poplist,
        patch_artist=True,
        medianprops={"color": "k", "linewidth": lw},
        whiskerprops={"color": "k"},
        capprops={"color": "k"},
        showfliers=False,
        flierprops={"c": "k", "markersize": 2})
    ax.set_ylabel(r'NodeAge', rotation=90, fontsize=12)
    ax.set_xticklabels(poplist, rotation=45, fontsize=8)
    # set list of colors
    colornames = list(mpl.colors.cnames.keys())[:len(poplist)]
    for patch, color in zip(box['boxes'], colornames):
        patch.set_facecolor(color)
        patch.set_linewidth(lw)
    fig.savefig("NodeDepth.pdf", bbox_inches="tight")
    return(None)


def getDivergence(infile, topos, pairs, minFreq, toposplot, pairsplot):
    """parses out distance file to return a distribution for each species pair
    """
    div_dict = defaultdict(lambda: defaultdict(list))
    tree_count = 0
    with open(infile,'r') as f:
        header = next(f).split()
        for line in f:
            tree_count += 1
            for i, div in enumerate(line.split()):
                topo, ind1, ind2 = header[i].split("_")
                pairname = "{}-{}".format(ind1, ind2)
                if div != "nan":  # should remove long list of nan causing issues with boxplot
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


def sumBranchLengths(infile, topos_subset, minFreq, nodedepthplot, step=10, topos=105):
    """
    """
    
    blen_box = []
    tree_count = 0
    topodict = {}
    with open(infile, 'r') as f:
        for line in f:
            if line.startswith("#"):
                x = line.split()
                topo = x[0][1:]
                leafs = [l.split("--")[1] for l in x[3:]]
                topodict[topo] = leafs
            else:
                tree_count += 1
                i = 0
                j = step
                blen = list(map(float, line.split()))
                blen_list = []
                while j <= len(blen):
                    blen_vals = blen[i+2:j]  # first 2 are outgroup root lengths
                    if np.isnan(blen_vals[0]):
                        blen_list.append(np.nan)
                    else:
                        topo_key = i//step
                        blen_topo = topodict["topo{}".format(topo_key+1)]
                        leaf_ix = [blen_topo.index(bt) for bt in blen_topo if bt.count('_') == 0]  # should return 3 or more
                        blen_max = max([blen_vals[bl] for bl in leaf_ix])  # get all values for leaf index and max
                        blen_maxix = blen_topo[blen_vals.index(blen_max)]  #  here is the leaf name 'Fun'
                        leaf_dist = [blen_vals[blen_topo.index(x1)] for x1 in blen_topo if blen_maxix in x1]  # now return all vals with instances of 'Fun' in topo
                        leaf_sum = sum(leaf_dist)
                        blen_list.append(leaf_sum)
                    i += step
                    j += step
                blen_box.append(blen_list)
    
    # boxplot of node depth
    data = list(zip(*blen_box))
    
    if nodedepthplot:
        boxplotB(data, topos_subset)
    
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
    if args.dist:
        getDivergence(infile, topos, pairs, args.minFreq, toposplot=args.tplot, pairsplot=args.pplot)
    elif args.blen:
        sumBranchLengths(infile, topos, args.minFreq, nodedepthplot=args.bplot)
