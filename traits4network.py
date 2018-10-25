#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 19:04:26 2018

@author: scott
"""
import sys

## for label
#sp_count = 0
#species = "aconitus"
#with open("CO1.Key.srt.txt", 'r') as co1:
#    for line in co1:
#        x = line.split()
#        if x[2] == species:
#            sp_count += 1
#            f.write("{}\t{}\t{}_{}\n".format(x[0], x[1], x[2], sp_count))
#        else:
#            sp_count = 1
#            f.write("{}\t{}\t{}_{}\n".format(x[0], x[1], x[2], 1))
#            species = x[2]
#
## prct ID cutoff
#with open("CO1.Key.srt.txt" ,'r') as co1:
#    for line in co1:
#        x = line.split()
#        if float(x[2]) > 94.00:
#            x[3] = x[1]
#            f.write("{}\n".format("\t".join(x)))
#        else:
#            x = line.split()
#            f.write("{}\n".format("\t".join(x)))
#
## replace fasta header
#with open("CO1.Key.srt.txt", 'r') as co1:
#    for line in co1:
#        x = line.split()
#        fastadict[x[0]] = x[-1]
#
#with open("CO1_Karama_28SEP18.clean.fa", 'r') as co1:
#    for line in co1:
#        if line.startswith(">"):
#            x = line.strip()[1:]
#            f.write(">{}\n".format(fastadict[x]))
#        else:
#            f.write(line)
#
## add lengths to key
#lendict = {}
#with open("CO1_Karama_28SEP18.clean.lab.fa.fai", 'r') as fai:
#    for line in fai:
#        x = line.split()
#        length = x[1]
#        lendict[x[0]] = x[1]
#
#with open("CO1.Key.srt.txt",'r') as co1:
#    for line in co1:
#        x = line.split()
#        x.append(lendict[x[3]])
#        f.write("{}\n".format("\t".join(x)))
#
## add groups from abgd
#groupdict = {}
#with open("CO1_Karama_28SEP18.part.9.txt", 'r') as group:
#    for line in group:
#        x = line.split()
#        for i in x[4:]:
#            groupdict[i]=x[1]
#
#with open("CO1.Key.srt.txt",'r') as co1:
#    for line in co1:
#        x=line.split()
#        try:
#            x.append(groupdict[x[3]])
#            f.write("{}\n".format("\t".join(x)))
#        except KeyError:
#            x.append("NAN")
#            f.write("{}\n".format("\t".join(x)))

# for traits
header = ["aconitus", "barbirostris", "nitidus", "peditaeniatus", "maculatus", "tessellatus", "culicifacies", "vagus","unknown"]
f = open("t", 'w')
with open(sys.argv[1], 'r') as log:
    log.next()  # skip header
    for line in log:
        try:
            sp = line.split()[0]
            headerlist = [0, 0, 0, 0, 0, 0, 0, 0, 0]
            headerlist[header.index(sp.split("_")[0])] += 1
            line = log.next()
            while line.startswith("\t"):
                spix = header.index(line.split()[0].split("_")[0])
                headerlist[spix] += 1
                line = log.next()
            f.write("{},{}\n".format(sp, ",".join(map(str, headerlist))))
        except StopIteration:
            break
f.close()
