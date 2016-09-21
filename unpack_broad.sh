#!/bin/bash
module load java/1.8
module load bio/picardtools/1.119

#for d in ./An*/; do
#    cd "$d"
#    for a in $(ls -1 *.tar.gz); do
#	tar -zxvf $a
#	done
#    rm *.gz
#    for f in ./*/;do
#	cd "$f"
#	bam=$(ls *.bam)
#	fastq_name=$(echo "$f" | cut -c 3- | sed 's/.$//')
#	SamToFastq I=${bam} F=${fastq_name}.1.fastq.tmp F2=${fastq_name}.2.fastq.tmp 
#	mv *.fastq.tmp ../
#	cd ..
#	done
#    cat *1.fastq.tmp > ${d}.1.fastq
#    cat *2.fastq.tmp > ${d}.2.fastq
#    paste <(paste - - - - <${d}.1.fastq) <(paste - - - - <${d}.2.fastq) | tr '\t' '\n' > ${d}.int.fastq
#    rm *.fastq.tmp
#    pigz -p 4 *.fastq
#    cd ..
#done
#for d in ./An*/; do
#    cd "$d"
    for f in ./*/;do
	cd "$f"
	bam=$(ls *.bam)
	fastq_name=$(echo "$f" | cut -c 3- | sed 's/.$//')
	SamToFastq I=${bam} F=${fastq_name}.1.fastq.tmp F2=${fastq_name}.2.fastq.tmp 
	mv *.fastq.tmp ../
	cd ..
    done
#    cat *1.fastq.tmp > ${d}.1.fastq
#    cat *2.fastq.tmp > ${d}.2.fastq
#    paste <(paste - - - - <${d}.1.fastq) <(paste - - - - <${d}.2.fastq) | tr '\t' '\n' > ${d}.int.fastq
#    #rm *.fastq.tmp
#    /afs/crc.nd.edu/user/s/ssmall2/bin/pigz -p 4 *.fastq
#    cd ..
#done