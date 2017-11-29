#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:32:39 2017
python liftover.py -v FOO.vcf -o FOO.vcf.new -t FOO.coordinate
1) run RepeatMasker on 2 genomes
2) align the 2 genomes using mugsy
3) convert species vcf to bed
3) MafFilter
    -remove full gaps
    -create liftover by passing bed to MafFilter
4) sort by chrom and coordinate sort -k5,5 -k7,7n
5) optional: include bedfile of all positions
     fai2bed.py
     bedtools getfasta -fi FASTA -bed pos.bed -tab -bedOut > FOO.bed
     # getfasta seems to return the wrong seq if it is not 50,80,100 line char
6) optional: convert freebayes format to gatk format using fb2gatk.py
print fx.__doc__
@author: stsmall
"""
from __future__ import print_function
import argparse
from collections import defaultdict

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
    with open(vcfRef, 'r') as vcf:
        for line in vcf:
            if line.startswith('##'):
                outStream.write(line)
                continue
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
        transfer.next()  # skip header
        for line in transfer:
            x = line.strip().split()
            chrom = x[0]
            orient = x[1]
            # pos_s = x[2]
            pos_e = x[3]
            newchrom = x[4]
            newpos_s = x[6]
            # newpos_e = x[7]
            transdict[chrom][pos_e] = (newchrom, newpos_s, orient)
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


def vcfformat(gt, formats, tri=False, invariant=False):
    """Reformat genotype field from 1 sample column in a vcf to orient with
    any repolarization of the reference and alternate alleles

    Parameters
    ------
    gt: list, list of genotype GT:AD:DP:GQ:PL 0/0:994,0:994:99:0,120,1800
    tri: bool, indicating that the site is tri-allelic
    invariant: bool, invariant genotypes are shorter

    Returns
    ------
    gt: list, modified list of genotype information
    """
    formats = formats.split(':')
    if invariant:
        try:
            if 'RGQ' in formats:
                # GT:AD:DP:RGQ from HaplotypeCaller
                ad = gt[1]
                dp = gt[2]
                gq = gt[3]
                pl = "500,500,0"
            else:
                # GT:AD:DP from UnifiedGenotyper
                ad = gt[1]
                dp = gt[2]
                gq = '99'
                pl = "500,500,0"
        except IndexError:
            import ipdb;ipdb.set_trace()
        gt = "{}:{}:{}:{}:{}".format(gt[0], ad, dp, gq, pl)
    elif tri:
        pass
        # TODO: triallelic
        # 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
        # possible error on incorrect number of PL sites
#        nalts = gt[1].count(",")
#        if nalts == 2:
#            # normal triallelic
#        elif nalts == 3:
#            # all 4 bases
        gt = ":".join(gt)
    else:
        if 'PGT' in formats or 'PID' in formats:
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
        else:
            # GT:AD:DP:GQ:PL
            ad = gt[1]
            dp = gt[2]
            gq = gt[3]
            pl = gt[-1]
            # reverse AD
            ad1, ad2 = ad.split(",")
            ad = "{},{}".format(ad2, ad1)
            # reverse PL
            pl1, pl2, pl3 = pl.split(",")
            pl = "{},{},{}".format(pl3, pl2, pl1)
    # rebuild genotype
        gt = "{}:{}:{}:{}:{}".format(gt[0], ad, dp, gq, pl)
    return(gt)


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
    # TODO: Condense to an algorithm rather than 'if' statements
    if "." in x[4]:
        # this site is invariant, so a fixed difference between genomes
        for i, sample in enumerate(x[9:]):
            gt = sample.split(":")
            gt[0] = gt[0].replace("0", "1")  # change all 0s to 1s
            geno = vcfformat(gt, x[8], invariant=True)  # change gt columns
            x[i + 9] = geno
        alt_a = x[3]
    elif (x[3] in alt_a) and (ref_a in x[4]):
        # the reference matches the alt allele and alt allele matches the ref
        # simply a case of derived being the majority allele in the reference
        if "," in x[4]:
            # triallelic, such a pain!
            aa = alt_a.split(",")
            r = aa.index(x[3])  # find the ref allele in the triallelic list
            if r == 0:
                # ref allele is 1st alt allele
                for i, sample in enumerate(x[9:]):
                    gt = sample.split(":")
                    if "0/0" in gt[0]:
                        gt[0] = "1/1"
                    elif "1/1" in gt[0]:
                        gt[0] = "0/0"
                    elif "2" in gt[0]:
                        if "1/2" in gt[0]:
                            gt[0] = "0/2"
                        elif "0/2" in gt[0]:
                            gt[0] = "1/2"
                    geno = vcfformat(gt, x[8], tri=True)
                    x[i + 9] = geno
            elif r == 1:
                # ref allele is 2nd alt allele
                for i, sample in enumerate(x[9:]):
                    gt = sample.split(":")
                    if "0/0" in gt[0]:
                        gt[0] = "2/2"
                    elif "2/2" in gt[0]:
                        gt[0] = "0/0"
                    elif "2" in gt[0]:
                        if "1/2" in gt[0]:
                            gt[0] = "0/1"
                        elif "0/2" in gt[0]:
                            gt[0] = "1/2"
                    geno = vcfformat(gt, x[8], tri=True)
                    x[i + 9] = geno
        else:
            # simply switch the genotypes to re-orient
            for i, sample in enumerate(x[9:]):
                gt = sample.split(":")
                if "0/0" in gt[0]:
                    gt[0] = '1/1'
                elif "1/1" in gt[0]:
                    gt[0] = '0/0'
                geno = vcfformat(gt, x[8])
                x[i + 9] = geno
    elif (x[3] in alt_a) and (ref_a not in x[4]):
        # the case where the ancestral allele does not appear in the vcf
        if "," in x[4]:
            # triallelic is a pain!
            for i, sample in enumerate(x[9:]):
                # change all 0 to 3
                gt = sample.split(":")
                gt[0] = gt[0].replace("0", "3")
                geno = vcfformat(gt, x[8], tri=True)
                x[i + 9] = geno
            alt_a = "{},{}".format(x[4], x[3])
        else:
            # build a triallelic site where there is no reference
            for i, sample in enumerate(x[9:]):
                # change all 0 to 1, all 1 to 2
                gt = sample.split(":")
                gt[0] = gt[0].replace("1", "2")
                gt[0] = gt[0].replace("0", "1")
                geno = vcfformat(gt, x[8])
                x[i + 9] = geno
            alt_a = "{},{}".format(x[3], x[4])
    elif (x[3] not in alt_a) and (ref_a in x[4]):
        if "," in x[4]:
            # triallelic is a pain!
            aa = x[4].split(",")
            r = aa.index(ref_a)  # find the ref allele in the triallelic list
            if r == 0:
                # change 1 to 0
                # change 0 to 1
                for i, sample in enumerate(x[9:]):
                    gt = sample.split(":")
                    if '1/1' in gt[0]:
                        gt[0] = '0/0'
                    elif '0/0' in gt[0]:
                        gt[0] = '1/1'
                    elif '2' in gt[0]:
                        if '0/2' in gt[0]:
                            gt[0] = '1/2'
                        elif '1/2' in gt[0]:
                            gt[0] = '0/2'
                    geno = vcfformat(gt, x[8], tri=True)
                    x[i + 9] = geno
                alt_a = "{},{}".format(x[3], x[4][1])
            elif r == 1:
                # change 2 to 0
                # change 0 to 2
                for i, sample in enumerate(x[9:]):
                    gt = sample.split(":")
                    if '2/2' in gt[0]:
                        gt[0] = '0/0'
                    elif '0/0' in gt[0]:
                        gt[0] = '2/2'
                    elif '1' in gt[0]:
                        if '0/1' in gt[0]:
                            gt[0] = '1/2'
                        elif '1/2' in gt[0]:
                            gt[0] = '0/1'
                    geno = vcfformat(gt, x[8], tri=True)
                    x[i + 9] = geno
                alt_a = "{},{}".format(x[4][0], x[3])
        else:
            for i, sample in enumerate(x[9:]):
                # change all 1s to 0s, all 0s to 1s
                gt = sample.split(":")
                if '0/0' in gt[0]:
                    gt[0] = '1/1'
                elif '1/1' in gt[0]:
                    gt[0] = '0/0'
                geno = vcfformat(gt, x[8])
                x[i + 9] = geno
            alt_a = x[3]
    else:
        if alt_a in x[4]:
            # change 0 to 2s
            for i, sample in enumerate(x[9:]):
                gt = sample.split(":")
                gt[0] = gt[0].replace("0", "2")
                geno = vcfformat(gt, x[8], tri=True)
                x[i + 9] = geno
            alt_a = "{},{}".format(x[4], x[3])
        elif (x[3] not in ref_a) and (x[3] not in alt_a) and (ref_a not in x[4]):
            for i, sample in enumerate(x[9:]):
                gt = sample.split(":")
                gt[0] = gt[0].replace("1", "2")
                gt[0] = gt[0].replace("0", "1")
                geno = vcfformat(gt, x[8], tri=True)
                x[i + 9] = geno
            alt_a = "{},{}".format(x[4], x[3])
        else:
            print("{}\t{}\t{}\n".format(ref_a, alt_a, x))
    return(x, alt_a)


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
    miss = 0
    print("executing liftover ...")
    tx = open("UnalignedCarryOver.bed", 'w')
    # t = open("nonmatchingref.out", 'w')
    with open(vcfFile, 'r') as vcf:
        for line in vcf:
            if not line.startswith("#"):
                x = line.strip().split()
                chrom = x[0]
                pos = x[1]
                try:
                    newchrom, newpos, orient = transdict[chrom][pos]
                    ref_a, alt_a = refdict[newchrom][newpos]
                    if orient == "-":
                        # print("Before: {} -- {}".format(ra, ref_a))
                        x[3] = reverseComplement(x[3])
                        x[4] = reverseComplement(x[4])
                        # print("After: {} -- {}".format(ra, ref_a))
                    if x[3] != ref_a:
                        # t.write("BEFORE\t{}\t{}\t{}\n".format(ref_a, alt_a, x))
                        miss += 1
                        x, alt_a = reorientGT(x, ref_a, alt_a)
                        # t.write("AFTER\t{}\t{}\t{}\n".format(ref_a, alt_a, x))
                    elif "." in x[4]:
                        for i, sample in enumerate(x[9:]):
                            gt = sample.split(":")
                            if len(gt) < 4:
                                # GT:AD:DP from UnifiedGenotyper
                                ad = gt[1]
                                dp = gt[2]
                                gq = '99'
                                pl = "0,500,500"
                            else:
                                # GT:AD:DP:RGQ from HaplotypeCaller
                                ad = gt[1]
                                dp = gt[2]
                                gq = gt[3]
                                pl = "0,500,500"
                            geno = "{}:{}:{}:{}:{}".format(gt[0], ad, dp, gq,
                                                           pl)
                            x[i + 9] = geno
                    x[1] = newpos
                    x[3] = ref_a
                    x[4] = alt_a
                    x[8] = "GT:AD:DP:GQ:PL"
                    outStream.write("{}\n".format("\t".join(x)))
                except KeyError:
                    tx.write("{}\t{}\n".format(x[0], x[1]))
    # t.close()
    tx.close()
    print("Mismatch Reference Allele {} times".format(miss))
    return(outStream)


if __name__ == "__main__":
    outstream = open(args.outFile, "w")
    outstream = headerVCF(args.vcfRef, outstream)
    refdict = loadVCF(args.vcfRef, args.refBed)
    transdict = loadTransfer(args.transfersFile)
    outstream = liftover(args.vcfFile, transdict, refdict, outstream)
    outstream.close()
