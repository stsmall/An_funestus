#!/bin/bash
#1 first pair
#2 second pair
#3 int outfile

gunzip -c $3 | paste - - - - - - - - | tee >(cut -f 1-4 | tr '\t' '\n' | gzip > $1) | cut -f 5-8 | tr '\t' '\n' | gzip -c > $2