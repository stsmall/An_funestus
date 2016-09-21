# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 13:32:07 2016

@author: stsmall
"""

from collections import defaultdict
trees=defaultdict(list)
with open("top.out",'r') as top:
    for line in top:
        if line.startswith("("):      
            try:
                while line.startswith("("):
                    trees[line.strip()].append(next(top))
                    line = next(top)
            except StopIteration:            
                pass

with open("mono.parvanlong.out",'r') as mono:
    with open("mono.tre",'r') as order:
        for line in order:
            trees[line.strip()].append(next(mono))
            
            
import re, subprocess         
start = 1
tree_line = 1
p=open("file",'w')
with open("file",'r') as trees:  #tree\tnumber\n
    for line in trees:
        x = line.strip().split()
        end = int(x[1])
        f=open("tmp.tre",'w')
        f.write(x[0])
        f.close()
        phylo_it = p.readlines()
        for phylo in phylo_it[start:end]:
            phylo.split(".")[2] #label
            command = "/afs/crc.nd.edu/user/s/ssmall2/programs_that_work/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -f J -m GTRGAMMA -p 12345 -s phy_reduced/aln." + phylo.split(".")[2] + ".phy.reduced -t tmp.tre -n boot." + phylo.split(".")[2] + " -o Anrivulorum_Kwa3inds"           
            proc = subprocess.Popen(command, shell=True)
            proc.wait()
            with open("RAxML_fastTreeSH_Support" + phylo.split(".")[2],'r') as part:
                f=open(str(tree_line)+".out",'w')                
                for tre in part:                
                    dp = re.findall(r'\[\d+\]', tre)
                    f.write(x[0]+"\naln." + phylo.split(".")[2] + ".phy.reduced\t %i\n" %dp)
                f.close()
        tree_line += 1
        start += end