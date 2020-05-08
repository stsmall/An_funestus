Snakefile
	1) make masks for data and training 
		use python makeMaskArmFilesFromQualAndRM.py from github.com/kern-lab/FILET/tree/master/maskingScripts
		required: gVCF from 2 pops, repeatGFF, reference fasta
		returns: fasta with masked ('N') sites
	2) make_filet_mask.py
	3) data2classify
		see examples at github.com/kern-lab/FILET/tree/master/maskingScripts
		required: phased fasta for each pop // arm (masked), anc fasta, bed file of steps
		returns: obs stat vectors for classification
Training sims
	abc_sims.py
Snakefile
	1) filet_statsMP.py, requires masks and associated C programs
	2) split : test/training sets
	3) training
	4) classify test sets
	5) classify data vectors

Previous efforts for this project:
################make mask file :: bedtools, makeFILETmask.py
# make mask file, see masking.txt
cat Van.neg-mask.srt.bed Fun.neg-mask.srt.bed | sort -k1,1 -k2,2n > FV.neg-mask.srt.bed

mergeBed -i FV.neg-mask.srt.bed > FV.neg-mask.srt.merge.bed

for chr in 2L 2R 3L 3R X;do                                                                                                                                                                                    
maskFastaFromBed -fi ${chr}.FV.fa -bed FV.neg-mask.srt.merge.bed -fo ${chr}.FV.mask.fa
done

# check options for script
for chr in 2L 2R 3L 3R X;do                                                                                                                                                                                    
python ssmall2/An_funestus/makeFILETmask.py -f ${chr}.FV.mask.fa -w 10000 -n 100
done

# make certain that mask is same length as sims
awk '/\/\// { found++ } found < 100001 { print}' FV.mask.fa.mask2 > FV.mask
