#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 17:15:28 2017

@author: scott
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('mafFile', metavar="mafFile", type=str,
                    help='path to maf file')
args = parser.parse_args()


def fixmafgap(mafFile):
    """
    """
    c = open("{}.changed".format(mafFile), 'w')
    f = open("{}.out".format(mafFile), 'w')
    with open(mafFile, 'r') as maf:
        for line in maf:
            if line.startswith("#"):
                f.write(line)
                c.write(line)
            elif line.startswith("a"):
                p1 = maf.next()
                p1ref = p1.strip().split()
                p2 = maf.next()
                p2ref = p2.strip().split()
                # s AgamP4.chr2R 8961 298 + 61545105 CAA
                if p2ref[3] == 0 or p1ref[3] == 0:
                    pass
                    c.write(line)
                    c.write(p1)
                    c.write(p2)
                elif p1ref[3] != p2ref[3]:
                    p1ref[6] = p1ref[6][0:-1]
                    p2ref[6] = p2ref[6][0:-1]
                    if int(p1ref[3]) > int(p2ref[3]):
                        p1ref[3] = str(int(p1ref[3]) - 1)
                    else:
                        p2ref[3] = str(int(p2ref[3]) - 1)
                    f.write(line)
                    f.write("{}\n".format(" ".join(p1ref)))
                    f.write("{}\n".format(" ".join(p2ref)))
                    c.write(line)
                    c.write("{}\n".format(" ".join(p1ref)))
                    c.write("{}\n".format(" ".join(p2ref)))
                else:
                    f.write(line)
                    f.write(p1)
                    f.write(p2)
            else:
                f.write(line)
    f.close()
    c.close()
    return(None)


if __name__ == '__main__':
    fixmafgap(args.mafFile)
