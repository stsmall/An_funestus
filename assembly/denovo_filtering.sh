#!/bin/bash
#denovo_filter fasta_base

#$1 name in current directory
#$2 path to discovar

#remove less than 1000bp
seqtk seq -L 1000 ${2}/a.lines.fasta > ${1}.1000.fasta

#make bed from faidx
samtools faidx ${1}.1000.fasta
awk '{if($2 >= 1000) print $1 "\t0\t" $2 "\t" $1}' ${1}.1000.fasta.fai > ${1}.selectSeq.1k.bed

#these should be the same
echo "total records less than 1000bp"
wc -l ${1}.1000.fasta.fai
wc -l ${1}.selectSeq.1k.bed

#megablastn, through qsub
blastn -db /scratch365/ssmall2/blast_db/db/nucleotide_nt/nt -num_threads 20 -outfmt '6 qseqid sscinames sseqid evalue bitscore qlen length pident mismatch qcovs qcovus qstart qend sstart send staxID' -max_target_seqs 1 -max_hsps 100 -query ${1}.1000.fasta > ${1}.discovar.megablastn.out

#sort and filter by evalue <.01 and take highest bit score
awk -F"\t" '$4 <= 0.01' ${1}.discovar.megablastn.out | sort -u -k1,1 | awk -F"\t" '$5 >= 50' > ${1}.discovar.megablastn.u.out

#histogram of non anopheline megablast
grep -v "Anopheles" ${1}.discovar.megablastn.u.out | cut -f2 | cut -d" " -f1 | sort | uniq -c | sort -n -k1,1 > ${1}.nonanopheles.genus.hist

#known anopheles
fgrep -w -e 'Anopheles' -e 'Aedes' -e 'Culex' -e 'Drosophila' ${1}.discovar.megablastn.u.out | grep -vw "Wolbachia" | cut -f1 | sort -u -k1,1 > ${1}.anopheles.mbn.contigs
fgrep -w -f ${1}.anopheles.mbn.contigs ${1}.selectSeq.1k.bed > ${1}.anopheles.bed

#known non-anopheles
fgrep -wv -f ${1}.anopheles.mbn.contigs ${1}.discovar.megablastn.u.out | cut -f 1 | sort -u -k1,1 > ${1}.nonanopheles.contigs

#unknown
cat ${1}.nonanopheles.contigs ${1}.anopheles.mbn.contigs | sort -u -k1,1 > ${1}.known.contigs
fgrep -wv -f ${1}.known.contigs ${1}.selectSeq.1k.bed > ${1}.unknown.bed

#echo "number of contigs in megablast"
sort -u -k1,1 ${1}.discovar.megablastn.out | wc -l
echo "not filtered out"
wc -l ${1}.discovar.megablastn.u.out
wc -l ${1}.anopheles.mbn.contigs
wc -l ${1}.nonanopheles.contigs

fastaFromBed -fi ${1}.1000.fasta -bed ${1}.unknown.bed -fo ${1}.unknown.fasta

#nucmer funestus
echo "nucmer funestus"
nucmer /afs/crc.nd.edu/user/s/ssmall2/ssmall2/references/Anfunestus.pacbio_merged.scf.fasta ${1}.unknown.fasta
rename out ${1}.unknown out.delta
delta-filter -q ${1}.unknown.delta > ${1}.q.delta
show-coords -THrlc ${1}.q.delta | awk '{if ($7 > 70 && $5 > 500) print $12"\t"$1"\t"$2"\t"$13"\t"$11}' > ${1}.nucmer.bed
cut -f4 ${1}.nucmer.bed | cut -d":" -f1 | sort -u -k1,1 > mapping.contigs
rm out.mgaps out.ntref

fgrep -w -f mapping.contigs ${1}.unknown.bed >> ${1}.anopheles.bed
fgrep -wv -f mapping.contigs ${1}.unknown.bed > ${1}.notmapping.bed
fastaFromBed -fi ${1}.1000.fasta -bed ${1}.notmapping.bed -fo ${1}.notmapping.fasta

wc -l ${1}.unknown.bed
wc -l mapping.contigs
wc -l ${1}.notmapping.bed
grep -c ">" ${1}.notmapping.fasta

#nucmer minimus
echo "nucmer minimus"
nucmer /afs/crc.nd.edu/user/s/ssmall2/ssmall2/references/anop_mini_minimus1.0.scaffolds.fasta ${1}.notmapping.fasta
rename out ${1}.unknown out.delta
delta-filter -q ${1}.unknown.delta > ${1}.q.delta
show-coords -THrlc ${1}.q.delta | awk '{if ($7 > 70 && $5 > 500) print $12"\t"$1"\t"$2"\t"$13"\t"$11}' > ${1}.nucmer.bed
cut -f4 ${1}.nucmer.bed | cut -d":" -f1 | sort -u -k1,1 >> mapping.contigs
sort -u -k1,1 mapping.contigs > ${1}.mapping.contigs
fgrep -w -f ${1}.mapping.contigs ${1}.unknown.bed >> ${1}.anopheles.bed
fgrep -wv -f ${1}.mapping.contigs ${1}.unknown.bed > ${1}.notmapping2.bed
fastaFromBed -fi ${1}.1000.fasta -bed ${1}.notmapping2.bed -fo ${1}.notmapping2.fasta

echo "results of nucmer"
wc -l ${1}.unknown.bed
wc -l ${1}.mapping.contigs
wc -l ${1}.notmapping2.bed
grep -c ">" ${1}.notmapping2.fasta
rm mapping.contigs

#blastn contigs that failed everything else; qsub
blastn -db /scratch365/ssmall2/blast_db/db/anopheles/anopheles.fsa -num_threads 20 -outfmt '6 qseqid sseqid evalue bitscore qlen length pident mismatch qcovs qcovus qstart qend sstart send' -max_target_seqs 1 -max_hsps 100 -query ${1}.notmapping2.fasta > ${1}.discovar.blastn.out

awk -F"\t" '$3 <= 0.01' ${1}.discovar.blastn.out | sort -u -k1,1 | awk -F"\t" '$4 >= 50' > ${1}.discovar.blastn.u.out
cut -f1 ${1}.discovar.blastn.u.out | cut -d":" -f1 > ${1}.anopheles.bn.contigs
fgrep -w -f ${1}.anopheles.bn.contigs ${1}.notmapping2.bed >> ${1}.anopheles.bed
fgrep -wv -f ${1}.anopheles.bn.contigs ${1}.notmapping2.bed > ${1}.nohits.bed

echo "results of blastn"
sort -u -k1,1 ${1}.discovar.blastn.out | wc -l
echo "not filtered out"
wc -l ${1}.discovar.blastn.u.out
wc -l ${1}.anopheles.bn.contigs
echo "contigs that are still unidentified"
wc -l ${1}.nohits.bed

sort -u -k1,1 ${1}.anopheles.bed > ${1}.anopheles.u.bed
wc -l ${1}.anopheles.u.bed
#redundans on all known anopheles contigs
fastaFromBed -fi ${1}.1000.fasta -bed ${1}.anopheles.u.bed -fo ${1}.redundans.in.fasta
samtools faidx ${1}.redundans.in.fasta
~/programs_that_work/redundans/redundans.py -v -i *.fastq.gz -f ${1}.redundans.in.fasta -t 12 -o ${1}.redundans --sspacebin /afs/crc.nd.edu/user/s/ssmall2/programs_that_work/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl --noscaffolding --nogapclosing