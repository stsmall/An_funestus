#!/bin/bash
#1 first pair
#2 second pair
#3 out.int
paste <(gunzip -c $1 | paste - - - -) <(gunzip -c $2 | paste - - - -) | tr '\t' '\n' | gzip -c > $3