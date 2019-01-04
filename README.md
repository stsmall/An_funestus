
# FASTA processing
## accessibility.py
build accessibility file similar to Malaria gen, 1000 genomes project. Requires pysam

## fai2bed.py
creates a bed file from a fai file. samtools faidx FOO.fa

## fasta2vcf.py
creates a fasta file from a single sample vcf

## fastarenamefromfilename.py
renames the header in the fasta file from the file name

## fillMAF.py
fills Ns in MAF file ...

## firstFasta.py
prints the first entry from a multifasta

## fixmafgap.py
fixes gaps in MAF file

## reverse_fasta.py
reverse complement a fasta file

## size_selection_fasta.sh
select only contigs of min len from assembly in fasta format

# Phylogenetic Trees

## calcD2stat.py
calculates the D2 and D1 statistics from Hibbins and Hahn 2018

## calcNodeHeights.py
summary of node heights. Requires ete3

## calcNodeHeights.py
summary of node heights. Requires ete3

## calcT1T2.allel.py
calculates the 1st and 2nd node times from Fontaine 2015. Requires scikit-allel

## maf2phy.awk
reformats MAF to phylip

## maf2raxml.py
reformats MAF to run in raxml

## makeUltrametric.py
reformats tree to be ultrameteric

## parseDistmat.py
parse a distmatrix for plotting

## pruneTips.py
prunes tips so only 1 per species. Requires ete3

## traits4network.py
formats a traits file for popart

## trees2reduce.py
a bunch of functions for phylogenetic trees

## tree_summ.py
calculates branch lengths


# BPP

## BppScrape.py
summarizes MAP from BPP output for plotting in R

## clustPhy.py
clusters phylip formatted files where each file contains 1 locus.

## makeBPPFile.py
creates input files for use with BPP

## makeCTL4BPP.py
make all control file for parallel run of BPP

## makeNegBed.py
creates an inverse bedfile. Used to select intergenic and intron using a gff


# ASSEMBLY
## assemblathon_stats.pl
script from Keith Bradnam to calculate summary stats for a fasta file

## deinterleave_gz.sh
fast deinterleave a read files

## denovo_filtering.sh
filtering on de novo assemblies following methods in Small et al. 2018

## interleave_gz.sh
opposite of deinterleave

## mauve2consensus.py
creates a consensus of multiple mauve alignments for gap filling

## sga-ice.notstupid.sh
changes sga-ice to not store tmp files

## w2rap.kmers.sh
iterative run of w2rap for different kmers

# VCF manipulations

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

# FILET

## makeFILETmask.py
creates a mask file for use with FILET when using masked data

## filet_clustering.py
clusters regions classified by FILET using a rolling mean

## ms2filet.py
use ms to create a a feature vector for use with FILET

# Other programs
## fastafillgeno.py
makes a complete geno format file for use with simonhmartin script. This is needed if you did not have invariant sites in the vcf.

## convertLD2cM.py
converts a file with Rho values to cM

## divTime.py
uses algorithm in https://doi.org/10.1101/281881 to calculate div time.

## estsfs_format.py
format data file for use in est-sfs from a vcf

## derivedVCF.py
build a polarized VCF file from the output of est-sfs

## crispr.py
example only

## format4phylonet.py
creates an input file for use with PhyloNet

## format_3PCLR.py
creates input file for use with 3P-CLR

## fish_circos.py
merges labels on a nucmer coords output and creates a circos formatted file

## popstats.py
file from pontussk to calculate population genetic statistics

## stairs2ms.py
converts stairwayPlot output to ms command

## Testimation2.nb
mathematica code for ttps://doi.org/10.1101/281881
