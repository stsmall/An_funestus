## PSL2tln.py
converts a PSL file from HAL tools to liftover.table file

## accessibility.py
creates a list of accessible sites following suggestions in Ag1000g by leveraging files produced by bedtools genomeCov and bam files

## assemblathon_stats.pl
calculate basic set of metrics from a genome assembly fasta. Author: Keith Bradnam, Genome Center, UC Davis

## deinterleave_gz.sh
fast deinterleave fastq files

## denovo_filtering.sh
filtering on de novo assemblies using blast and nucmer against known species genomes

## fish_circos.py
creates a circos formatted file from nucmer output of genome alignments

## interleave_gz.sh
opposite of deinterleave on fastq files

## liftover.py
reorients a vcf file to coordinates of an aligned reference

## mauve2consensus.py
I made this script to create a consensus between 2 de novo genome assemblies
each from a difference individual but from the same species. The de novo
assemblies were then scaffolded in ragout. N's in 1 scaffolded assembly were
somtimes present as sequence in the other. If you align the 2 genomes with
Mauve to product a xmfa, you can then produce a consensus with fewer gaps than
either alone

## outgroupConsensus.py
creates a single consensus genotype in a VCF from multiple outgroups. usage to detect polarize conflicts among outgroups and when programs only allow for a single
sample to be designated as outgroup.

## remove_dups.tln.py
removes duplicates, paralogs, from liftover.table produced by hal

## sga-ice.notstupid.sh
Iterative error correction with sga-ice version which does not store tmp files and fills your drive

## size_select_fasta.py
simple stupid script to size select on a fasta file using an index. It is hard coded so probably just use seqtk instead.

## w2rap.kmers.sh
iterative run of w2rap for testing different kmers
