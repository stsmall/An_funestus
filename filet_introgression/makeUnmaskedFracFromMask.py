#!/usr/bin/env python
import sys

maskFile = sys.argv[1]
f = open(f"{maskFile}-unmaskedFrac", 'w')
with open(maskFile, 'r') as m:
    for line in m:
        if line.strip():
            if not line.startswith("//"):
                frac = 0
                while line.strip():
                    x = line.split()
                    frac += float(x[2]) - float(x[1])
                    line = next(m)
                f.write("{}\n".format(1-frac))
f.close()
