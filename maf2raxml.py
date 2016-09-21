# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:32:40 2016
maf2raxml.py maffile -s SIZE -f FASTA_list -e exclude -p partition_size

this should take a maf file and select for the minimum size alignment
remove any individual rows: maftools
stitch to a multifasta: maftools
fasta2phy: bash w/ RaxML
create a partition for RaxMl
run RaxML
concatenate trees
@author: stsmall
"""
#requires: RaxML (threaded version), script fasta2phy.sh (included with RaxML under usefulScripts), maf2fasta.pl from mugsy install

import argparse
import os
import subprocess
import glob

def get_args():
    parser = argparse.ArgumentParser(description='Parses MAF file from mugsy, runs RaxML')  
    parser.add_argument('-s','--align_size', help='minimum alignment size to save',type=int,default=10000)      
    parser.add_argument('-p','--partition_size', help='partition size for RaxML',type=int,default=10000)   
    parser.add_argument('-m','--maf_file',required=True, help='location of maf file') 
    parser.add_argument('-t','--threads',help='number of threads for RaxML',type=int,default=10)   
    parser.add_argument('-rx','--raxml',help='path to RaxML')
    parser.add_argument('-mf2fa','--maf2fasta',help='path to mugsy maf2fasta')    
    parser.add_argument('-swd','--workingdir',help='path to working directory',default=os.getcwd())    
    args = parser.parse_args()
    return args

def select_maf(maf,align_size):                  
    with open("trunc_maf.maf",'w') as f:
        with open(maf,'r') as maf_in:
            for line in maf_in:
                if line.startswith("a"):
                    mult = line.split()[3].split("=")[1]
                    a = line.split()
                    if int(mult) >= 11:
                        line = next(maf_in)                        
                        if line.startswith("s"):
                            block_len = line.split()[3]
                            if int(block_len) >= align_size:
                                f.write("%s\n" %" ".join(a))
                                try:
                                    while line.startswith("s"):
                                        f.write(line)
                                        line = next(maf_in)
                                except StopIteration:            
                                    pass
                                
                f.write("\n")
                             
                
def convertMaf2phy():
    command = "grep label trunc_maf.maf > label.list"
    print command            
    proc = subprocess.Popen(command, shell=True)
    proc.wait()
    with open("label.list",'r') as lab:
        for line in lab:
            label = line.split()[2].split("=")[1]
            command = "perl maf2fasta.pl " + label + " < trunc_maf.maf > " + "aln."+label+".fa" 
            print command            
            proc = subprocess.Popen(command, shell=True)
            proc.wait()
    for aln in glob.glob("*.fa"):
        label = aln.split(".")[1]        
        command = "bash convertFasta2Phylip.sh " + aln + " > aln." + label + ".phy" 
        print command            
        proc = subprocess.Popen(command, shell=True)
        proc.wait()

def phy2part(part_size):
    for phy in glob.glob("*.phy"):    
        label = phy.split(".")[1]
        with open("aln.{}.part".format(label),'w') as f:
            with open(phy,'r') as fl:
                first_line = fl.readline()
                part_end = first_line.split()[1]
                i = 1
                j = int(part_size)
                count = 1 
                while j < int(part_end):
                   f.write("DNA,p%i=%i-%i\n" %(count,i,j)) 
                   count += 1
                   i = i + int(part_size)
                   j = j + int(part_size)
                   if j >= int(part_end):
                       count += 1
                       j = int(part_end)
                       f.write(("DNA,p%i=%i-%i\n" %(count,i,j)))
                       
def runraxml(raxml):
    for phy in glob.glob("*.phy"):
        label = phy.split(".")[1]        
        command = raxml + "/raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -p 12345 -f K -s " + "aln." + label + ".phy -n rax." + label + ".tre -T 10 -o Anrivulorum_Kwa3inds -q aln." + label + ".part -M -t starting.tre"                 
        print command            
        proc = subprocess.Popen(command, shell=True)
        proc.wait()
        
def sumtrees(summary_topo,part_tre):        #really this should sort by unique topology, currently uses non-unique whcih makes more trees
    '''sorts all trees into the summary topologies, takes a summary file that is constructed
        by removing all branch length information with sed -e "s/\:[0-9]*.[0-9]*//g" Anfun.part.tre | sort -u and then using sort -nu'''
    from collections import defaultdict
    import re
    topo_store=defaultdict(list)
    with open(summary_topo,'r') as topo:
        for line in topo:
            topo_store[line.strip()]
    with open(part_tre,'r') as tre:
        for line in tre:
            no_bl = re.sub(r':\d*\.\d{3,}','',line.strip())
            #log_bl = re.sub(r':\d*\.\d{3,}','',line.strip()) 
            topo_store[no_bl].append(line.strip())
    cnt = 1
    for topo in topo_store.keys():
        f = open("top.{}.tre".format(cnt),'w')
        for tree in topo_store[topo]:        
            f.write(tree+"\n")
        f.close()
        cnt += 1
        
def collapsetrees(string):
    '''I loaded a set of non-unique trees these into densitree then exported the trees, it made a nexus tree file that contained only unique trees,
    I then renamed the numbers back to the newick format, although there is prob a converter somewhere. Then I ran ete3 compare against the list and recorded
    the number of times I saw a 0 robinson score, I then grouped by the monophyletic funestus proper'''
    #treedist in phanghorn in R
    #phybin will provide a matrix of distance, where 0 is the same tree
    #ete3 compare -t -r --taboutput #look for score of 0 is same tree
    #dendropy, also has sumtrees
    """Generate parenthesized contents in string as pairs (level, contents)."""
    stack = []
    for i, c in enumerate(string):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            start = stack.pop()
            yield (len(stack), string[start + 1: i])      
def main(): 
    args = get_args()
    os.chdir(args.workingdir)
    select_maf(args.maf_file,args.align_size)
    convertMaf2phy(args.maf2fasta,args.raxml)
    phy2part(args.partition_size)
    runraxml(args.raxml)
    #all in 1 file for densitree    
    command = "cat *bestTree* > all.tre"
    print command            
    proc = subprocess.Popen(command, shell=True)
    proc.wait()
    #all topologies in 1 file for sort, uniq -c    
    command = "sed -e 's/:0.[0-9]*//g' all.tre > topo.tre"
    print command            
    proc = subprocess.Popen(command, shell=True)
    proc.wait()
    
if __name__ == '__main__':
    main()