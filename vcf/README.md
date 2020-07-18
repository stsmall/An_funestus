## addAncAllele2Ref.py
Create an ancestral fasta file from a reference by replacing the ref base 
with the ancestral base. Requires Biopython. 
Ancestral file format has 3 columns
CHROM POS ALLELE

## chunkvcf.py
breaks vcf into smaller chunks for faster processing
`chunkvcf.py --vcf FOO.vcf --integer 10000`

## derivedVCF.py
creates a VCF with polarization defined by est-sfs (Keightly and Jackson 2018)

## estsfs_format.py
Create input files for est-sfs (Keightly and Jackson 2018)
file format should be from vcftools --counts
`estsfs_format.py -i ingroup.counts -o1 outgroup.counts [-o2 outgroup.counts -o3 outgroup.counts]`

## fb2gatk.py
convert freebayes vcf format to that used by gatk

## fixfbmiss.py
fixes freebayes missing to be in the same format as gatk

## fixgtfieldlength.py
The fields can be different lengths if the site is physically phased in gatk

## fixmissing.py
reformats missing sites in VCf

## gatknorm.py
normalizes a gatk file for merging

## geno2vcf.py
creates a geno file (Plink) from a vcf

## normgatkfields.py
normalizes gatk fior merging

## vcf2bedgraph.py
create a bedgraph from a vcf file

## vcf2fasta_consensus.py

## vcf2linkSelThin.py
thin a vcf file to reduce the effects of linked selection.

## vcf2sample_resort.py
remake the header for a vcf

## vcf_fill_gvcf.py
fill a gvcf file with 0 band
