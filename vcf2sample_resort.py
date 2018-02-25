#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 14:04:05 2017
@author: scott

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, {}}'
$46 ,$19 ,$162 ,$139 ,$122 ,$158 ,$43 ,$188 ,$63 ,$166 ,$160 ,$153 ,$142 ,$171 ,$79 ,$128 ,$138 ,$108 ,$131 ,$123 ,$41 ,$98 ,$99 ,$145 ,$102 ,$65 ,$14 ,$184 ,$146 ,$38 ,$75 ,$163 ,$36 ,$183 ,$125 ,$55 ,$124 ,$13 ,$47 ,$148 ,$193 ,$190 ,$176 ,$34 ,$157 ,$89 ,$32 ,$101 ,$152 ,$112 ,$92 ,$156 ,$111 ,$12 ,$87 ,$62 ,$27 ,$72 ,$109 ,$107 ,$80 ,$150 ,$71 ,$51 ,$42 ,$11 ,$25 ,$180 ,$137 ,$129 ,$29 ,$120 ,$84 ,$168 ,$169 ,$69 ,$82 ,$24 ,$15 ,$33 ,$81 ,$78 ,$172 ,$28 ,$164 ,$141 ,$10 ,$35 ,$60 ,$96 ,$100 ,$22 ,$91 ,$30 ,$134 ,$74 ,$113 ,$178 ,$186 ,$151 ,$76 ,$103 ,$192 ,$119 ,$50 ,$23 ,$49 ,$167 ,$73 ,$40 ,$64 ,$136 ,$144 ,$59 ,$149 ,$165 ,$48 ,$155 ,$177 ,$195 ,$17 ,$161 ,$21 ,$37 ,$116 ,$18 ,$95 ,$26 ,$140 ,$110 ,$194 ,$106 ,$147 ,$52 ,$70 ,$68 ,$31 ,$154 ,$104 ,$191 ,$53 ,$133 ,$77 ,$126 ,$130 ,$45 ,$66 ,$189 ,$58 ,$135 ,$83 ,$61 ,$115 ,$39 ,$86 ,$57 ,$173 ,$9 ,$94 ,$90 ,$182 ,$159 ,$174 ,$44 ,$85 ,$20 ,$181 ,$16 ,$179 ,$143 ,$118 ,$132 ,$105 ,$114 ,$187 ,$67 ,$185 ,$117 ,$97 ,$175 ,$88 ,$170 ,$121 ,$93 ,$56 ,$54 ,$127
"""
import sys

with open(sys.argv[1], 'r') as vcf:
    for line in vcf:
        if line.startswith("#CHROM"):
            x = line.strip().split()
            samples = x[9:]
            newindex = sorted(samples)
            samplelist = []
            for s in newindex:
                samplelist.append(" ${}".format(1+x.index(s)))
            print("{}".format(",".join(samplelist)))
            break