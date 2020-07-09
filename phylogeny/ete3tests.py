#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 15:25:10 2018

@author: scott
"""


from ete3 import PhyloTree
import numpy as np

outgroup = "rivulorum_F790"
t = PhyloTree("((((rivulorum_F790:0.25862,((vaneedeni_KwaF773:0.0,vaneedeni_KwaF774:0.0,vaneedeni_KwaF775:0.0,vaneedeni_KwaF780:0.0,vaneedeni_KwaF782:0.0,vaneedeni_KwaF783:0.0,vaneedeni_KwaF786:0.0):0.00055,vaneedeni_KwaF784:0.00055)0.982:0.02123)0.176:0.00899,((funestus_TanF561:0.00055,((funestus_MozF804:0.00055,funestus_Ugf401:0.00055)0.000:0.00055,(funestus_GhaF265:0.00055,(funestus_GhaF264:0.0,funestus_Ken4590:0.0,funestus_MozF123:0.0,funestus_MozF260:0.0,funestus_MozF29:0.0,funestus_MozF35:0.0,funestus_TanF601:0.0,funestus_Ugf399:0.0,funestus_Ugf403:0.0,funestus_Zam281:0.0):0.00055)0.000:0.00055)0.856:0.00120)0.891:0.00241,(funestuscf_MALAF105_7:0.00055,(funestuscf_MALAF99_4:0.0,funestuscf_MALF98_2:0.0):0.00055)0.993:0.00860)0.854:0.00228)0.966:0.00732,((parensis_KwaF761:0.0,parensis_KwaF762:0.0,parensis_KwaF766:0.0,parensis_KwaF767:0.0,parensis_KwaF768:0.0,parensis_KwaF769:0.0,parensis_KwaF835:0.0):0.00053,parensis_KwaF851:0.00055)0.982:0.00089)0.948:0.00384,((longipalpusC_11:0.0,longipalpusC_13:0.0):0.00055,longipalpusC_551_12533:0.00076)1.000:0.00051,((longipalpusC_15:0.0,longipalpusC_16:0.0,longipalpusC_551_12634:0.0):0.00055,(longipalpusC_4:0.00055,longipalpusC_12:0.00055)1.000:0.00055)0.709:0.00055);")
t.set_species_naming_function(lambda node: node.name.split("_")[0])  # n.species, n.name
t.set_outgroup( t&outgroup )
#taxon1 = ["parensis", "longipalpusC", "vaneedeni"]
#taxon2 = ["funestus", "funestuscf", "vaneedeni"]
taxon = [["parensis", "longipalpusC", "vaneedeni"], ["funestus", "funestuscf", "vaneedeni"]]
taxdict = {}
for i, tax in enumerate(taxon):
    nodesupport = []
    nodeage = []
    for t in treelist:
        if t.check_monophyly(values=tax, target_attr="species"):
            samples = []
            for sp in tax:
                samples.extend(t.search_nodes(species=sp))
            ancnode = t.get_common_ancestor(samples)
            nodeage.append(ancnode.dist)
            nodesupport.append(ancnode.support)
    taxdict[i] = (nodeage, nodesupport)



    if not winlist:
        winarray = np.ones(len(treelist), dtype=bool)
    mtreelist, winarray = getMonophyletic(treelist, quart, winarray)
    btreelist, winarray = supportFilt(mtreelist, quart, winarray)
    if nodes:
        nh1, nh2 = nodeHeights(btreelist, quart)
    else:
        nh1 = []
        nh2 = []
    return(treelist, winarray, nh1, nh2)




def nodeHeights(mtreelist, quart):
    """Calculates the node heights in a set of trees.

    Parameters
    ------
    trees: ete3 object, returned from function loadtrees
    quart: list, list of topologies to calculate node heights
    windows: list, genome coordinates for which the trees were made

    Returns
    ------
    nh1: float, average value of T1 for topology
    nh2: float, average value of T2 for topology
    """
    print("calculating node heights for quartets: {}".format(quart))
    nh1 = []
    nh2 = []
    for t in mtreelist:
        P1 = t.search_nodes(species=quart[0])
        P2 = t.search_nodes(species=quart[1])
        P3 = t.search_nodes(species=quart[2])
        for p1, p2 in combinations(P1, P2, 2):
            nh2.append(t.get_distance(p1, p2))
        for p1, p3 in combinations(P1, P3, 2):
            nh1.append(t.get_distance(p1, p3))
        for p2, p3 in combinations(P2, P3, 2):
            nh1.append(t.get_distance(p2, p3))
    print("T1: {}, T2: {}".format(np.mean(np.array(nh1)),
                                  np.mean(np.array(nh2))))
    return(nh1, nh2)


def getMonophyletic(treelist, quart, winarray):
    """
    """
    p1, p2, p3, p4 = quart
    mtreelist = []
    for i, t in enumerate(treelist):
        if cMono(t, p1) and cMono(t, p2) and cMono(t, p3) and cMono(t, p4):
            Out = t.get_common_ancestor(t.search_nodes(species=quart[-1]))
            t.set_outgroup(Out)
            mtreelist.append(t)
        else:
            winarray[i] = False
#            t2 = t.collapse_lineage_specific_expansions()
#            print(t2)
#            for node in t.split_by_dupes():
#                print(node)
#            ntrees, ndups, sptrees = t.get_speciation_trees()
#            for spt in sptrees:
#                print(spt)
    return(mtreelist, winarray)