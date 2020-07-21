# Coalescent times
## calcT1T2.allel.py
calculates the 1st and 2nd node times from Fontaine 2015 from a tree with mulitple individuals per species. Requires scikit-allel

## T1T2Stats.py
Calculates the coalescent times of a triplet following the formula in Fontaine et al 2015 using NewickTools.

## calcD1D2stat.py
calculates the D2 and D1 statistics from Hibbins and Hahn 2018. Can subsamble from a phylogeny with multiple indiviuals per species and calculate standard errors by iterating along the genome in windows.

## newickStats.py
rewrite of D1D2 statistics from Hibbins and Hahn 2018 with a corrected AND logic, also some comment notes. Written for direct comparison to the calcD1D2stat.py script

# Tree Stats
## calcTreeStats.py
Calculate average support, root ages, monophyly, rogue taxa

## pruneTips.py
prunes tips from a phylogenetic tree to only 1 individual (haplotype) per species. Requires ete3

## tree_summ.py
calculates branch lengths

## trees2RFdistance.py
create pairwise matrix of robison-foulds distance between phylogenetic trees. requires ete3

## calcRootDivAstral.py
Returns the time to the MRCA of a ASTRAL tree.

## makeUltrametric.py
reformats tree to be ultrameteric. requires ete3

# Local MSC trees
## astral_pp.py
astral_pp.py is a program to construct local species trees from collections
of gene trees distributed along the genome. This idea is similar to that
of Thawornwattana et al 2018, where they used BPP for this same thing. ASTRAL
is much faster and therefore it is easier to explore the clustering number for
fine scale analysis.

## gff2fastaAln.py
This module creates input files for the methods of Thawornattana 2018.
Using a gff file it parses out individual alignments of exons and non-coding
sequences.
`python gff2fastaAln.py --gff FOO.gff --aln FOO.aln.fa --chromlen int [--minCDS] [--prct] [--distance int] [--mxlength int] [--mnlength] [--bpp] [--clust int]`

## makeCTL4BPP.py
Takes the output of makeBPPFile.py and creates control (input) files for running BPP in parallel

## bpp_scrape.py
Programs to collect data from mulitple run of BPP and format into a weights file for chromosome painting.
`bpp_scrape.py -p nCDS.bpp. -s .out -c 100000 --scafs 2L 2R 3L 3R X
returns: Coordinate File, Twisst-style weights file`

## clustPhy.py
clusters phylip formatted files into a single file with many loci. Input phylip are assumed to have only 1 locus.

## create_loci_cluster.py
This code is designed to take a set of fasta or phylip files each with 1 locus
and cluster them into 1 file containing many loci. This specific iteration of
the code was designed to generate input files for BPP. The file name is expected to carry information on the loci $chrom.$cds_type.$coords.aln.33.fa = 2L.cds.102218-102369.aln.33.fa
`python create_loci_clusters.py --aln_path tests/fasta_files --aln_format fa`

## makeBPPFile.py
Superseded by gff2fastaAln.py
`python makeBPPFile.py --gff gff --fasta FOO.fa [--exons] [--distance int] [--length int] [--chromlen int]`

# Parsing
## parseDistmat.py
parse a distmatrix for plotting

## parseTwisst.py
will parse a divergence and branch length file produced by twisst https://github.com/simonhmartin/twisst
`parseTwisst.py -i dist.txt -t [topos] -p [pairs]`

## parse_twisst_dist_len.py
Finding the tree heights or coalescent times between species pairs is best done using ETE3. This is simple for trees with 1 allele per species but takes much more time for trees with >1 allele per species. The program twisst can estimate distances and branch lengths when calculating weights.This script will parse the Dist and Len files produced by twisst. Namely specifying topologies or frequencies. twisst https://github.com/simonhmartin/twisst
`twisst.py -t TREES_in -w WEIGHTS_out -D DISTANCE_out -L BRANCH_LEN_out --outputTopos TOPO_out`

# Format
## format4phylonet.py
Creates a newick file with command and headers for running PhyloNet.

## maf2phy.awk
reformats MAF to phylip. Author: Bernhard Haubold, haubold@evolbio.mpg.de

## maf2raxml.py
reformats MAF to run in raxml
`maf2raxml.py maffile -s SIZE -f FASTA_list -e exclude -p partition_size`

## rename.AFCmtDNA.msOut.py
Hard coded. Renames tips in ms-format tree with ACF species names

## traits4network.py
formats a traits file for POPART networks. some hard coding