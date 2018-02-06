#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 13:57:52 2017

@author: stsmall
"""
from __future__ import print_function
from __future__ import division
# from IPython.display import HTML
import numpy as np
from ete3 import PhyloTree
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', "--treefile", type=str, required=True,
                    help="treefile in newick, 1 per line")
parser.add_argument('-g', "--groups", nargs='+', action='append',
                    required=True, help="quartet of species to calculate,"
                    " assumes form: P1 P2 P3. can be given multiple times")
parser.add_argument('-w', "--windows", type=str,
                    help="coordinates for each tree")
parser.add_argument("--dlm", type=str, default="_",
                    help="delimeter denoting species")
parser.add_argument("-o", "--outgroup", type=str,
                    help="outgroup species for rooting")
args = parser.parse_args()


def LoadTrees(treefile, quart, outgroup, dlm):
    """Reads and stores phylogenetic trees from a file

    Parameters
    ------
    treefile: file, file of newick trees, 1 per line
    outgroup: str, last entry from quartet

    Returns
    ------
    treelist: obj, ete3 object of trees

    """
    print("loading trees...")
    treelist = []
    with open(treefile, 'r') as newick:
        for line in newick:
            if not line.startswith("NA"):
                t = PhyloTree(line)
                t.set_species_naming_function(lambda node: node.name.split(dlm)[0])
                t.set_outgroup( t&outgroup )
                treelist.append(t)
    return(treelist)


def WriteTrees(treelist):
    """Rewrite reformatted or rerooted trees

    Parameters
    ------
    treelist: list, list of tree objects

    Returns
    ------
    file: file, writes out tree objects in nexus format

    """
    f = open("rooted_trees.nex", 'w')
    for t in treelist:
        f.write("{}\n".format(t.write()))
    f.close()
    return(None)


def cMono(tree, taxon):
    """Checks if samples are monophyletic in tree
    """
    return(tree.check_monophyly(values=taxon, target_attr="species")[0])


def rogueTree(t, tax):
    """
    """
    import ipdb; ipdb.set_trace()
    v = t.get_monophyletic(values=tax, target_attr="species")
    m = list(v)
    if len(m) <= 2:
        s = t.get_common_ancestor(m)
        tips = s.get_leaf_names()  # list of all individuals, find the odd one
        uniqtips = np.unique(tips, return_counts=True)
        species = uniqtips[0]
        inds = uniqtips[1]

        # print tree, individual name
    else:
        pass
    # tree.remove_child(child)
    # tree.prune(nodes, preserve_branch_length=True)
    # t2 = t.collapse_lineage_specific_expansions()
    return(None)


def AgeAndSupport(treelist, taxon):
    """Calculates the support and node age if groups in taxon are monophyletic
    """
    taxdict = {}
    agelist = []
    supportlist = []
    for tax in taxon:
        nodesupport = []
        nodeage = []
        for t in treelist:
            if cMono(t, tax):
                v = t.get_monophyletic(values=tax, target_attr="species")
                ancnode = list(v)[0]
                nodeage.append(ancnode.get_farthest_leaf()[1])
                nodesupport.append(ancnode.support)
#                nodeage.append(np.mean([t.get_distance(ancnode, i) for i in t.iter_descendants()]))
            else:
                nodeage.append(np.nan)
                nodesupport.append(np.nan)
                # rogueTree(t, tax)
        agelist.append(nodeage)
        supportlist.append(nodesupport)
    taxdict["age"] = agelist
    taxdict["support"] = supportlist
    return(taxdict)


def FilterTree(treelist):
    """Return the MRCA of all leafs/tips
    """
    sum_support = []
    root_height = []
    for t in treelist:
        root_height.append(t.get_farthest_leaf()[1])  # filter short trees
        supp = []
        for node in t.get_descendants():
            if len(node) > 1:
                supp.append(node.support)
        sum_support.append(np.mean(supp))
    return(root_height, sum_support)


def AgeStats(taxdict, quart):
    """Plot histogram and stats for branch lengths and support
    """
    # age hist
    for pop, age in enumerate(taxdict['age']):
        ma = np.nanmean(age)  # mean age
        nt = np.count_nonzero(~np.isnan(age))  # num trees
        print("\n{0}, numtrees {1}, meanage {2:.6f}".format(quart[pop], nt, ma))
        label = "".join([j[0] for j in quart[pop]])
        age = np.array(age)
        age = age[~np.isnan(age)]
        plt.hist(age, bins=30, alpha=0.5, label=label)
    plt.legend(loc='upper right')
    plt.ylabel("node age")
    plt.savefig("NodeHeight.pdf")
    plt.clf()
    return(None)


def SupportStats(taxdict, quart):
    # support hist
    for pop, supp in enumerate(taxdict['support']):
        ms = np.nanmean(supp)  # mean age
        nt = np.count_nonzero(~np.isnan(supp))  # num trees
        print("\n{0}, numtrees {1}, meansupport {2:.2f}".format(quart[pop], nt, ms))
        label = "".join([j[0] for j in quart[pop]])
        supp = np.array(supp)
        supp = supp[~np.isnan(supp)]
        plt.hist(supp, bins=30, alpha=0.5, label=label)
    plt.legend(loc='upper right')
    plt.ylabel("support")
    plt.savefig("NodeSupport.pdf")
    plt.clf()
    return(None)


def WindowStats(windows, taxdict, quart, root_height, sum_support):
    """
    """
    winlist = []
    f = open("window_stats.tsv", 'w')
    r = open("tree_stats.tsv", 'w')
    with open(windows, 'r') as win:
        header = win.next()
        for line in win:
            winlist.append(line.rstrip("\n"))
    label = []
    for tax in quart:
        label.append("".join([j[0] for j in tax]))
    f.write("Summary\t{}\t{}\n".format(header.rstrip("\n"), "\t".join(label)))
    r.write("{}\trootHeight\tmeanSupport\n".format(header.rstrip("\n")))
    nage = zip(*taxdict['age'])
    nsupp = zip(*taxdict['support'])
    try:
        for i, age in enumerate(nage):
            f.write("NodeHeight\t{}\t{}\n".format(winlist[i], "\t".join(map(str, age))))
            f.write("NodeSupport\t{}\t{}\n".format(winlist[i], "\t".join(map(str, nsupp[i]))))
            r.write("{}\t{}\t{}\n".format(winlist[i], root_height[i], sum_support[i]))
    except TypeError:
        import ipdb;ipdb.set_trace()
    f.close()
    r.close()
    return(None)


if __name__ == "__main__":
    quart = args.groups
    treelist = LoadTrees(args.treefile, quart, args.outgroup, args.dlm)
    WriteTrees(treelist)
    taxdict = AgeAndSupport(treelist, quart)
    AgeStats(taxdict, quart)
    SupportStats(taxdict, quart)
    root_height, sum_support = FilterTree(treelist)
    if args.windows:
        WindowStats(args.windows, taxdict, quart, root_height, sum_support)

##test tree
#t='(rivulorum_F790:0.1188,(((longipalpusC_551_12634:6e-09,(longipalpusC_13:6e-09,(longipalpusC_16:6e-09,longipalpusC_551_12533:6e-09)0.921:6e-09)0.909:5e-09)0.014:0.00222105,(longipalpusC_15:6e-09,(longipalpusC_12:6e-09,(longipalpusC_11:0.00112065,(longipalpusC_4:6e-09,(parensis_KwaF762:5e-09,((parensis_KwaF761:6e-09,parensis_KwaF766:6e-09)0.92:6e-09,(parensis_KwaF767:0,parensis_KwaF768:0,parensis_KwaF769:0,parensis_KwaF835:0,parensis_KwaF851:0)1:6e-09)0:6e-09)0.583:6e-09)0.292:0.00110003)0:2.27e-07)0:5e-09)0.711:2.305e-05)0.955:0.00882347,(((((vaneedeni_KwaF782:6e-09,vaneedeni_KwaF780:6e-09)0.921:6e-09,vaneedeni_KwaF774:6e-09)0:6e-09,(vaneedeni_KwaF784:6e-09,vaneedeni_KwaF783:6e-09)0.767:6e-09)0.889:5e-09,(vaneedeni_KwaF775:6e-09,(vaneedeni_KwaF773:6e-09,vaneedeni_KwaF786:0.00112541)0.367:5e-09)1:7.08e-07)0.995:0.0102093,((funestuscf_MALAF105_7:0,funestuscf_MALAF99_4:0,funestuscf_MALF98_2:0)1:0.00447164,((funestus_MozF123:6e-09,((((funestus_MozF35:0,funestus_MozF804:0,funestus_Zam281:0)1:6e-09,funestus_TanF561:6e-09)0.936:5e-09,funestus_TanF601:6e-09)0.395:6e-09,funestus_MozF29:6e-09)0.405:0.00334382)0.646:6e-09,(funestus_GhaF264:6e-09,(funestus_Ken4590:6e-09,(funestus_GhaF265:6e-09,(funestus_Ugf399:6e-09,(funestus_Ugf403:6e-09,(funestus_MozF260:6e-09,funestus_Ugf401:6e-09)0.731:5e-09)0.85:6e-09)0.133:0.00222583)0.459:6e-09)0:5e-09)0.453:6e-09)0.894:0.00222803)0.789:6e-09)0.726:0.00356974)1:0.1188);'
#tree = PhyloTree(t)
#tree.set_species_naming_function(lambda node: node.name.split("_")[0])
#tree.set_outgroup( tree&'rivulorum_F790')
#
#tree.check_monophyly(["longipalpusC"], target_attr="species")
## 0 is bool
## 2 is problem nodes
#tree.get_monophyletic(values=["longipalpusC"], target_attr="species")
#tree.remove_child(child)
#tree.prune(nodes, preserve_branch_length=True)

