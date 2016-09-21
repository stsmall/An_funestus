#!/bin/bash
set -e
set -o pipefail

#seg 's/Anparensis_Kwa8ind/Anrivulorum_Kwa3inds/g' sga-ice.notstupid.sh

### Create 'ec' directory for the final error-corrected files ###
mkdir -p /scratch365/ssmall2/Anparensis_Kwa8inds/ec
cd /scratch365/ssmall2/Anparensis_Kwa8inds/ec
#pigz -p 40 -d /scratch365/ssmall2/Anparensis_Kwa8inds.nocontaminants.{1,2}.fastq.gz

### Link input files to temp dir ###
ln -s /scratch365/ssmall2/Anparensis_Kwa8inds/Anparensis_Kwa8inds.nocontaminants.{1,2}.fastq .

### Run preprocessing for each input file ###
echo '### Start sga preprocessing ###'
sga preprocess --no-primer-check --pe-mode 1 --permute-ambiguous --min-length 0 --out Anparensis_Kwa8inds.nocontaminants.pp.fastq Anparensis_Kwa8inds.nocontaminants.1.fastq Anparensis_Kwa8inds.nocontaminants.2.fastq
echo '    Build new index from the preprocessed input files'
sga index -a ropebwt -t 40 -p all.ec.k0 *.pp.fastq
for f in *pp*; do mv $f $(echo $f | sed 's/.pp./.ec.k0./g') ; done
echo '### Finished sga preprocessing ###'

### Do k-mer based correction
echo '### Starting 8 rounds of k-mer based error correction with the following ks: 40, 60, 100, 125, 150, 200, 225, 250  ###'
echo '  Running correction round with k=40 '
echo '    Correcting file Anparensis_Kwa8inds.nocontaminants.ec.k0.fastq'
sga correct --count-offset 2 --threads 40 --prefix all.ec.k0 --outfile Anparensis_Kwa8inds.nocontaminants.ec.k40.fastq -k 40 --learn Anparensis_Kwa8inds.nocontaminants.ec.k0.fastq
rm *ec.k0.*

echo '    Build new index'
sga index -a ropebwt -t 40 -p all.ec.k40 *ec.k40.fastq
echo '  Running correction round with k=60 '
echo '    Correcting file Anparensis_Kwa8inds.nocontaminants.ec.k40.fastq'
sga correct --count-offset 2 --threads 40 --prefix all.ec.k40 --outfile Anparensis_Kwa8inds.nocontaminants.ec.k60.fastq -k 60 --learn Anparensis_Kwa8inds.nocontaminants.ec.k40.fastq
rm *ec.k40.*

echo '    Build new index'
sga index -a ropebwt -t 40 -p all.ec.k60 *ec.k60.fastq
echo '  Running correction round with k=100 '
echo '    Correcting file Anparensis_Kwa8inds.nocontaminants.ec.k60.fastq'
sga correct --count-offset 2 --threads 40 --prefix all.ec.k60 --outfile Anparensis_Kwa8inds.nocontaminants.ec.k100.fastq -k 100 --learn Anparensis_Kwa8inds.nocontaminants.ec.k60.fastq
rm *ec.k60.*

echo '    Build new index'
sga index -a ropebwt -t 40 -p all.ec.k100 *ec.k100.fastq
echo '  Running correction round with k=125 '
echo '    Correcting file Anparensis_Kwa8inds.nocontaminants.1.ec.k100.fastq'
sga correct --count-offset 2 --threads 40 --prefix all.ec.k100 --outfile Anparensis_Kwa8inds.nocontaminants.ec.k125.fastq -k 125 --learn Anparensis_Kwa8inds.nocontaminants.ec.k100.fastq
rm *ec.k100.*

echo '    Build new index'
sga index -a ropebwt -t 40 -p all.ec.k125 *ec.k125.fastq
echo '  Running correction round with k=150 '
echo '    Correcting file Anparensis_Kwa8inds.nocontaminants.ec.k125.fastq'
sga correct --count-offset 2 --threads 40 --prefix all.ec.k125 --outfile Anparensis_Kwa8inds.nocontaminants.ec.k150.fastq -k 150 --learn Anparensis_Kwa8inds.nocontaminants.ec.k125.fastq
rm *ec.k125.*

echo '    Build new index'
sga index -a ropebwt -t 40 -p all.ec.k150 *ec.k150.fastq
echo '  Running correction round with k=200 '
echo '    Correcting file Anparensis_Kwa8inds.nocontaminants.ec.k150.fastq'
sga correct --count-offset 2 --threads 40 --prefix all.ec.k150 --outfile Anparensis_Kwa8inds.nocontaminants.ec.k200.fastq -k 200 --learn Anparensis_Kwa8inds.nocontaminants.ec.k150.fastq
rm *ec.k150.*

echo '    Build new index'
sga index -a ropebwt -t 40 -p all.ec.k200 *ec.k200.fastq
echo '  Running correction round with k=225 '
echo '    Correcting file Anparensis_Kwa8inds.nocontaminants.ec.k200.fastq'
sga correct --count-offset 2 --threads 40 --prefix all.ec.k200 --outfile Anparensis_Kwa8inds.nocontaminants.ec.k225.fastq -k 225 --learn Anparensis_Kwa8inds.nocontaminants.ec.k200.fastq
rm *ec.k200.*

echo '    Build new index'
sga index -a ropebwt -t 40 -p all.ec.k225 *ec.k225.fastq
echo '  Running correction round with k=250 '
echo '    Correcting file Anparensis_Kwa8inds.nocontaminants.ec.k225.fastq'
sga correct --count-offset 2 --threads 40 --prefix all.ec.k225 --outfile Anparensis_Kwa8inds.nocontaminants.1.ec.k250.fastq -k 250 --learn Anparensis_Kwa8inds.nocontaminants.1.ec.k225.fastq
rm *ec.k225.*
echo '### Finished all k-mer based correction rounds ###'

### Do overlap-based correction
echo '### Run overlap-based correction ###'
echo '    Build new index from the output of k-mer correction'
sga index -a ropebwt -t 40 -p all.ec.k250 *ec.k250.fastq
echo '    Overlap-correction of file Anparensis_Kwa8inds.nocontaminants.1.ec.k250.fastq'
sga correct --prefix all.ec.k250 --algorithm overlap --error-rate 0.010000 -m 40 --threads 40 --outfile Anparensis_Kwa8inds.nocontaminants.int.final.ecOv.fastq Anparensis_Kwa8inds.nocontaminants.ec.k250.fastq
echo '### Finished overlap-based correction ###'

### After overlap correction, .fastq files are converted to .fasta because the quality values in the .fastq file produced by 'sga correct' not always match the length of the reads if insertions/deletions were corrected.
### Newer versions of sga may solve this issue
echo '### After overlap correction, convert fastq to fasta files ###'
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; print ">"header, seq}' < Anparensis_Kwa8inds.nocontaminants.int.final.ecOv.fastq > Anparensis_Kwa8inds.nocontaminants.int.final.ecOv.fastq.fasta

pigz -p40 *final.ecOv.fastq
pigz -p40 *ec.k250.fastq
echo; echo; echo '### ALL DONE.  Final error-corrected files after k-mer-based and overlap-based correction are /afs/crc.nd.edu/user/s/ssmall2/ssmall2/Anparensis_Kwa8inds/ec/*.final.ecOv.f*q* ###'