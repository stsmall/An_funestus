from Bio import SeqIO
from Bio import Seq
import re

pattern = [NNNNNGG, CCNNNNNNN]

pam = re.compile(pattern)

crisprdict = {}
fasta_sequences = SeqIO.parse(fastaFile, 'fasta')
for fasta in fasta_sequences:
	header, sequence = fasta.id, str(fasta.seq)
	# sequence.reverse_complement()
	crispr = [(m.span(), m.group()) for m in re.finditer(pam, sequence)]
    crisprdict[h] = crispr

f = open("crispr.out", 'w')
for chr in crisprdict.keys():
	for c in crisprdict[chr]:
		pos, seq = c
		f.write('{}\t{}\t{}\t{}\n'.format(chr, pos[0], pos[1], seq))
f.close()