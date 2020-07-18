## fai2bed.py
creates a bed file from a fai file. 
`samtools faidx FOO.fa`
`python fai2bed.py FOO.fa.fai  # > pos.bed`
`bedtools getfasta -fi FOO.fa -bed pos.bed -tab -bedOut > FOO.bed`

## fasta2phy.py
Convert fasta alignemnts to relaxed phylip ones in constant memory.
Written by Lucas Sinclair; James B. Pease
`fasta2phy.py IN.fasta OUT.phylip`

## fasta2raxml.py
Written by Jon Palmer (2016) nextgenusfs@gmail.com
Script runs Mafft -> trimAl -> RAxML. takes multi-fasta as input

## fastarenamefromfilename.py
renames the header in the fasta file from the file name

## firstFasta.py
prints the first entry from a multifasta

## gff2fastaAln.py
Using a gff file it parses out individual alignments of exons and non-coding
sequences
`python gff2fastaAln.py --gff FOO.gff --aln FOO.aln.fa [--exons] [--distance int] [--length int] [--chromlen int]`

## reverse_fasta.py
reverses sections of a fasta file given by coordinates. Useful for correcting small artifact inversions. Some hard-coding.

##sfsFromFasta.py
returns the site frequency spectrum from a fasta file
`sfsFromFasta.py --fasta FOO.fa --ancestral FOO.anc.fa`
