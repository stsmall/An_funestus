
fasta=$1 #basename with no .fa

samtools faidx ${fasta}.fasta
awk '{if($2 >= 500) print $1 "\t0\t" $2 "\t" $1}' ${fasta}.fasta.fai > ${fasta}.selectSeq500.bed
fastaFromBed -fi ${fasta}.fasta -bed ${fasta}.selectSeq500.bed -name -fo ${fasta}.500.fasta