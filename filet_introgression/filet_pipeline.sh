# run sims and generate stats
qsub
	model.run # sims
	model.proc.run # stats 

# make ROC and training sets
cd $sp
mkdir ROC
for ms in *fvout;do                                                                                                                                    
    sed -n '1,10001p' $ms > fvec
    sed -n -e '1p' -e '10001,$p' $ms > ${ms}.2
    paste ../2L.cols fvec > ROC/${ms}.ROC.fvec
done
mv mig12.msout3.fvout.2 mig12.msOut
mv mig21.msout3.fvout.2 mig21.msOut
mv noMig.msout3.fvout.2 noMig.msOut
python ~/programs_that_work/FILET/buildThreeClassTrainingSet.py trainingSims/$sp/ trainingSets/threeClass.${sp}.mask.fvec

# train classifier
qsub
	python ~/programs_that_work/FILET/trainClassifier.py trainingSets/threeClass.${sp}.mask.fvec all classifiers/threeClass.${sp}.mask.p

# determine cut off from ROC w/ classifier	
#python ~/programs_that_work/FILET/classifyChromosome.py classifier/threeClass.${sp}.mask.p featureVectors/slide/ .1 results/${sp}
