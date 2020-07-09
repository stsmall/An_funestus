## PSL2tln.py
converts a PSL file from HAL tools to liftover.table file

## remove_dups.tln.py
removes duplicates, paralogs, from liftover.table

## vcf2bedgraph.py
create a bedgraph from a vcf file

## vcf2linkSelThin.py
thin a vcf file to reduce the effects of linked selection.

## vcf2sample_resort.py
remake the header for a vcf

## vcf_fill_gvcf.py
fill a gvcf file with 0 band

## fixmissing.py
reformats missing sites in VCf

## addAncAllele2Ref.py
take bedfile of ancestral allele positions and mutates a fasta to ref allele. Requires Biopython

## chunkvcf.py
breaks vcf into chunks for faster processing

## fb2gatk.py
convert freebayes vcf format to that used by gatk

## fixfbmiss.py
fixes freebayes missing to be in the same format as gatk

## fixgtfieldlength.py
The fields can be different lengths if the site is physically phased in gatk

## gatknorm.py
normalizes a gatk file for merging

## geno2vcf.py
creates a geno file (Plink) from a vcf

## liftover.py
reorients a vcf file to coordinates of a reference

## normgatkfields.py
normalizes gatk fior merging

## outgroupConsensus.py
creates a single consensus genotype from multiple outgroups
