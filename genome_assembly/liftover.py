#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:32:39 2017
python liftover.py -v1 FOO.vcf -v2 REF.vcf -o NAME -t FOO.tln
1) run RepeatMasker on 2 genomes
2) align the 2 genomes using ProgCactus
3) convert species vcf to bed; vcf2bedgraph.py
3) MafFilter
    -remove full gaps; Alnfilter2
    -liftover; LiftOver
4) sort by chrom and coordinate sort -k5,5 -k7,7n
5) gatkNorm.py; to normalize format fields
7) run liftover.py
    * optional: include bedfile of all positions w/ -b option
         fai2bed.py
         bedtools getfasta -fi FASTA -bed pos.bed -tab -bedOut > FOO.bed
         # getfasta seems to return the wrong seq if it is not
         50,80,100 line char
    * optional: convert freebayes format to gatk format using fb2gatk.py
8) bcftools sort
9) bgzip & bcftools index -t
10) bcftools merge

*filter vcf for before, remove '*' and indels, vcfstreamsort before merge
print fx.__doc__
@author: stsmall
"""
from __future__ import print_function
import argparse
from collections import defaultdict
import ipdb
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-v1', "--vcfFile", type=str, action="store",
                    required=True, help="name of infile")
parser.add_argument('-v2', "--vcfRef", type=str, action="store", required=True,
                    help="name of reference vcf")
parser.add_argument('-b', "--refBed", type=str, action="store", required=False,
                    help="if all sites not in VCFref then pass bedfile")
parser.add_argument('-o', "--outFile", type=str, action="store",
                    help="name of outfile")
parser.add_argument("-t", "--transfersFile", action='store', required=True,
                    help="Chromosome and position transfer table")
args = parser.parse_args()


def headerVCF(vcfRef, outStream):
    """Write header from input vcf to output vcf

    Parameters
    ------
    vcfFile: file, file in vcf format
    outstream: file, writeable file for new vcf

    Returns
    ------
    outstream: file, open file with new header

    """
    print("writing header ...")
    if "gz" in vcfRef:
        with gzip.open(vcfRef, 'rb') as vcf:
            for line in vcf:
                if line.startswith(b'##'):
                    outStream.write(line)
                else:
                    break
    else:
        with open(vcfRef, 'r') as vcf:
            for line in vcf:
                if line.startswith('##'):
                    outStream.write(line)
                else:
                    break
    return(outstream)


def loadVCF(vcfRef, refBed):
    """Load ref and alternate alleles from vcfRef for verification

    Parameters
    -----
    vcfRef: file, vcf formatted file for the new coordinate system
    refbed: file, bed formatted file needed if all invariants are not in vcf

    Returns
    ------
    refdict: dict, default dictionary where keys are also dictionary. Format is
        'chrom' : 'pos' : ('ra', 'aa')

    """
    print("loading ref vcf ...")
    refdict = defaultdict(dict)
    if "gz" in vcfRef:
        with gzip.open(vcfRef, 'rb') as vcf:
            for line in vcf:
                if not line.startswith(b"#"):
                    x = line.split()
                    chrom = x[0]
                    pos = x[1]
                    ra = x[3]
                    aa = x[4]
                    refdict[chrom][pos] = (ra, aa)
    else:
        with open(vcfRef, 'r') as vcf:
            for line in vcf:
                if not line.startswith("#"):
                    x = line.split()
                    chrom = x[0]
                    pos = x[1]
                    ra = x[3]
                    aa = x[4]
                    refdict[chrom][pos] = (ra, aa)
    if refBed:
        if "gz" in refBed:
            with gzip.open(refBed, 'rb') as bed:
                for line in bed:
                    x = line.strip().split()
                    chrom = x[0]
                    pos = x[2]
                    ra = x[3]
                    try:
                        refdict[chrom][pos]
                    except KeyError:
                        refdict[chrom][pos] = (ra, ".")
        else:
            with open(refBed, 'r') as bed:
                for line in bed:
                    x = line.strip().split()
                    chrom = x[0]
                    pos = x[2]
                    ra = x[3]
                    try:
                        refdict[chrom][pos]
                    except KeyError:
                        refdict[chrom][pos] = (ra, ".")
    return(refdict)


def loadTransfer(transfersFile):
    """Load a file with translated coordinates.

    Parameters
    ------
    transferFile: file, Transfers file should have seven fields for oldChrom,
        orientation, oldStart, oldEnd, newChrom, orientation, newStart, newEnd.

    Returns
    ------
    transdict: dict, default dictionary where keys are also dictionary. Format
        is 'chrom' : 'pos' : ('newchrom', 'newpos', 'orientation')

    """
    print("loading transfer file ...")
    transdict = defaultdict(dict)
    with open(transfersFile, 'r') as transfer:
        for line in transfer:
            if "+" not in line.split()[1]:
                line = transfer.next()
            x = line.strip().split()
            chrom = x[0]
            orient_r = x[1]  # this should be positive for the ref
            if orient_r != "+":
                raise Exception("Reference/Query Orientation should be + orientation")
            pos_e = x[3]
            newchrom = x[4]
            orient_t = x[5]  # this can be positive or negative
            newpos_e = x[7]
            try:
                transdict[chrom][pos_e]  # does this exist, likely no
                print("Key Exists")
            except KeyError:
                transdict[chrom][pos_e] = (newchrom, newpos_e, orient_t)
    return(transdict)


def reverseComplement(base):
    """Return reverse complement of a nucleotide when the orientation is neg
    strand.

    Parameters
    ------
    base: str, one of nucleotide base A,T,G,C

    Returns
    ------
    base: str, reverse complement of nucleotide base A,T,C,G

    """
    if "," in base:
        # multiple alternative alleles
        base1 = base.split(",")
        for i, b in enumerate(base1):
            if b == 'A':
                base = 'T'
            elif b == 'T':
                base = 'A'
            elif b == 'C':
                base = 'G'
            elif b == 'G':
                base = 'C'
            base1[i] = base
        base = ",".join(base1)
    else:
        if base == 'A':
            base = 'T'
        elif base == 'T':
            base = 'A'
        elif base == 'C':
            base = 'G'
        elif base == 'G':
            base = 'C'
    return(base)


def reorientGT(x, ref_a, alt_a):
    """Repolarize GT columns based on matching the ref or alt sites if the ref
    allele between vcf and refvcf do not match. Any time there is a diff
    between the 2 reference genomes.

    Parameters
    ------
    x: list, list containing 1 line from a vcf
    ref_a: str, the reference allele from the reference vcf
    alt_a: str, the alternate allele if any, from the reference vcf

    Returns
    ------
    x: list, modified list of line from vcf
    alt_a: modified alternate allele for the new vcf

    """
    formats = x[8].split(":")
    try:
        if "." in x[4]:
            # this site is invariant, so a fixed difference between genomes
            for i, sample in enumerate(x[9:]):
                gt = sample.split(":")
                if ("./." in gt[0]) or (gt[0] == '.'):
                    x[i + 9] = "./.:.:.:.:."
                elif (".|." in gt[0]):
                    x[i + 9] = ".|.:.:.:.:."
                else:
                    gt[0] = gt[0].replace("0", "1")  # change all 0s to 1s
                    ad = gt[formats.index('AD')]
                    dp = gt[formats.index('DP')]
                    gq = gt[formats.index('GQ')]
                    pl = gt[formats.index('PL')]
                    geno = "{}:0,{}:{}:{}:500,500,{}".format(gt[0], ad, dp, gq, pl)
                    x[i + 9] = geno
            x[4] = x[3]
        elif x[3] in alt_a:
            if ref_a in x[4]:
                for i, sample in enumerate(x[9:]):
                    gt = sample.split(":")
                    if ("./." in gt[0]) or (gt[0] == '.'):
                        x[i + 9] = "./.:.:.:.:."
                    elif (".|." in gt[0]):
                        x[i + 9] = ".|.:.:.:.:."
                    else:
                        if "0/0" in gt[0] or "0|0" in gt[0]:
                            gt[0] = gt[0].replace("0", "1")
                        elif "1/1" in gt[0] or "1|1" in gt[0]:
                            gt[0] = gt[0].replace("1", "0")
                        # fix formating
                        ad = gt[formats.index('AD')]
                        dp = gt[formats.index('DP')]
                        gq = gt[formats.index('GQ')]
                        pl = gt[formats.index('PL')]
                        # reverse AD
                        ad1, ad2 = ad.split(",")
                        ad = "{},{}".format(ad2, ad1)
                        # reverse PL
                        pl1, pl2, pl3 = pl.split(",")
                        pl = "{},{},{}".format(pl3, pl2, pl1)
                        geno = "{}:{}:{}:{}:{}".format(gt[0], ad, dp, gq, pl)
                        x[i + 9] = geno
                x[4] = x[3]
                # x[3] = x[4]
            elif ref_a not in x[4]:
                # build a triallelic site where there is no reference
                for i, sample in enumerate(x[9:]):
                    # change all 0 to 1, all 1 to 2
                    gt = sample.split(":")
                    if ("./." in gt[0]) or (gt[0] == '.'):
                        x[i + 9] = "./.:.:.:.:."
                    elif (".|." in gt[0]):
                        x[i + 9] = ".|.:.:.:.:."
                    else:
                        gt[0] = gt[0].replace("1", "2")
                        gt[0] = gt[0].replace("0", "1")
                        ad = gt[formats.index('AD')]
                        dp = gt[formats.index('DP')]
                        gq = gt[formats.index('GQ')]
                        pl = gt[formats.index('PL')]
                        ad1, ad2 = ad.split(",")
                        ad = "0,{},{}".format(ad1, ad2)
                        pl1, pl2, pl3 = pl.split(",")
                        pl = "500,500,{},500,{},{}".format(pl1, pl2, pl3)
                        geno = "{}:{}:{}:{}:{}".format(gt[0], ad, dp, gq, pl)
                        x[i + 9] = geno
                x[4] = "{},{}".format(x[3], x[4])
            else:
                ipdb.set_trace()
        elif x[3] not in alt_a:
            # case 2
            if ref_a in x[4]:
                for i, sample in enumerate(x[9:]):
                    # change all 1s to 0s, all 0s to 1s
                    gt = sample.split(":")
                    if ("./." in gt[0]) or (gt[0] == '.'):
                        x[i + 9] = "./.:.:.:.:."
                    elif (".|." in gt[0]):
                        x[i + 9] = ".|.:.:.:.:."
                    else:
                        if '0/0' in gt[0]:
                            gt[0] = gt[0].replace("0", "1")
                        elif '1/1' in gt[0]:
                            gt[0] = gt[0].replace("1", "0")
                        ad = gt[formats.index('AD')]
                        dp = gt[formats.index('DP')]
                        gq = gt[formats.index('GQ')]
                        pl = gt[formats.index('PL')]
                        ad1, ad2 = ad.split(",")
                        ad = "{},{}".format(ad1, ad2)
                        pl1, pl2, pl3 = pl.split(",")
                        pl = "{},{},{}".format(pl3, pl2, pl1)
                        geno = "{}:{}:{}:{}:{}".format(gt[0], ad, dp, gq, pl)
                        x[i + 9] = geno
                x[4] = x[3]
            elif alt_a in x[4]:
                # build a triallelic site where there is no reference
                for i, sample in enumerate(x[9:]):
                    gt = sample.split(":")
                    if ("./." in gt[0]) or (gt[0] == '.'):
                        x[i + 9] = "./.:.:.:.:."
                    elif (".|." in gt[0]):
                        x[i + 9] = ".|.:.:.:.:."
                    else:
                        # change all 0 to 1, all 1 to 2
                        gt[0] = gt[0].replace("1", "2")
                        gt[0] = gt[0].replace("0", "1")
                        ad = gt[formats.index('AD')]
                        dp = gt[formats.index('DP')]
                        gq = gt[formats.index('GQ')]
                        pl = gt[formats.index('PL')]
                        ad1, ad2 = ad.split(",")
                        ad = "0,{},{}".format(ad1, ad2)
                        pl1, pl2, pl3 = pl.split(",")
                        pl = "500,500,{},500,{},{}".format(pl1, pl2, pl3)
                        geno = "{}:{}:{}:{}:{}".format(gt[0], ad, dp, gq, pl)
                        x[i + 9] = geno
                x[4] = "{},{}".format(x[3], x[4])
        elif "." in alt_a:
            for i, sample in enumerate(x[9:]):
                gt = sample.split(":")
                if ("./." in gt[0]) or (gt[0] == '.'):
                    x[i + 9] = "./.:.:.:.:."
                elif (".|." in gt[0]):
                    x[i + 9] = ".|.:.:.:.:."
                else:
                    # change all 0 to 1, all 1 to 2
                    gt[0] = gt[0].replace("1", "2")
                    gt[0] = gt[0].replace("0", "1")
                    ad = gt[formats.index('AD')]
                    dp = gt[formats.index('DP')]
                    gq = gt[formats.index('GQ')]
                    pl = gt[formats.index('PL')]
                    ad1, ad2 = ad.split(",")
                    ad = "0,{},{}".format(ad1, ad2)
                    pl1, pl2, pl3 = pl.split(",")
                    pl = "500,500,{},500,{},{}".format(pl1, pl2, pl3)
                    geno = "{}:{}:{}:{}:{}".format(gt[0], ad, dp, gq, pl)
                    x[i + 9] = geno
            x[4] = "{},{}".format(x[3], x[4])
        else:
            print("{}".format("\t".join(x)))
            print("{}\t{}".format(ref_a, alt_a))
    except ValueError:
        ipdb.set_trace()
    return(x)


def liftover(vcfFile, transdict, refdict, outStream):
    """Performs liftover between vcfs

    Parameters
    ------
    vcfFile: file, vcf to liftover
    transdict: dict, dictionary of transfer coords, returned from loadTransfer
    refdict: dict, dictionary from reference vcf, returned from loadVCF
    outStream: file, file to write new liftover vcf

    Returns
    ------
    outStream: file, completed file

    """
    refmatch = 0
    refmismatch = 0
    unaligned = 0
    reffix = 0
    triallelic = 0
    invariant = 0
    variant = 0
    print("executing liftover ...")
    tx = open("UnalignedCarryOver.bed", 'w')
    t = open("ChangedSites.out", 'w')
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                outStream.write(line)
            elif not line.startswith("#"):
                x = line.strip().split()
                chrom = x[0]
                pos = x[1]
                try:
                    newchrom, newpos, orient = transdict[chrom][pos]
                    ref_a, alt_a = refdict[newchrom][newpos]

                    # check if invariant in Ref
                    if alt_a == ".":
                        invariant += 1  # invariant in ref
                    else:
                        variant += 1  # variant in ref

                    # change allele orientation by strand
                    if orient == "-":
                        if x[4] == ".":
                            x[3] = reverseComplement(x[3])
                        else:
                            x[3] = reverseComplement(x[3])
                            x[4] = reverseComplement(x[4])

                    # Ref_refAllele and Lift_refAllele
                    # case where they do not match
                    if x[3] != ref_a:

                        refmismatch += 1
                        forline = "\t".join((x[0], x[1], x[3], x[4], x[9]))
                        t.write("Before\n{}\n".format(forline))

                        if "," in x[4]:
                            triallelic += 1
                            x[4] = 'NA'  # does not handle tri-allelic
                        else:
                            x = reorientGT(x, ref_a, alt_a)

                        forline = "\t".join((x[0], x[1], ref_a, x[4], x[9]))
                        t.write("After\n{}\n".format(forline))
                        x[0] = newchrom
                        x[1] = newpos
                        x[3] = ref_a
                        if "NA" not in x[4]:
                            reffix += 1
                    # case where they are the same
                    elif x[3] == ref_a:
                        refmatch += 1
                        x[0] = newchrom
                        x[1] = newpos
                    else:
                        raise Exception("RefAllele does not dif or equal liftoverRef Allele")
                    if 'NA' not in x[4]:
                        outStream.write("{}\n".format("\t".join(x)))
                except KeyError:
                    unaligned += 1
                    tx.write("{}\t{}\n".format(x[0], x[1]))  # unaligned
    tx.close()  # close unaligned
    t.close()  # close changed
    print("Matching reference:{}".format(refmatch))
    print("Mismatch reference:{}".format(refmismatch))
    print("Mismatch fixed:{}".format(reffix))
    print("Mismatch not fixed:{}".format(refmismatch - reffix))
    print("Unaligned:{}".format(unaligned))
    print("Fixed Differences:{}".format(invariant))
    print("Polylmorphic:{}".format(variant))
    print("Skipped {} Triallelic".format(triallelic))
    return(outStream)


if __name__ == "__main__":
    outstream = open(args.outFile, "w")
    outstream = headerVCF(args.vcfRef, outstream)
    refdict = loadVCF(args.vcfRef, args.refBed)
    transdict = loadTransfer(args.transfersFile)
    ipdb.set_trace()
    outstream = liftover(args.vcfFile, transdict, refdict, outstream)
    outstream.close()
