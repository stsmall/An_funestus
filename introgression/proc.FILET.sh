#!/bin/bash
# bash proc.FILET.sh mig12.sims.msOut FV.mask 50000 12 10 22 ~/programs_that_work/FILET 10000
# bash proc.FILET.sh sim mask reps pop1 pop2 path window mig {1,2,3}
sim=$1
mask=$2
loci=$3
p1=$4
p2=$5
p3=$6
FILET=$7
window=$8

grep -n -A1 "msmove" $sim | sed -n 's/^\([0-9]\{1,\}\).*/\1d/p' | sed -f - $sim > ${sim}2
head -n2 $sim > ${sim}.msh
sed -i "s/msmove $p3 1/msmove $p3 $loci/g" ${sim}.msh
cat ${sim}.msh ${sim}2 > ${sim}3
rm -f ${sim}.msh
rm -f ${sim}2
cat ${sim}3 | $FILET/msMaskAllRows $mask | python $FILET/removeNedOutColumnsFromMsFileKeepStars.py stdin > ${sim}3.mask
if [ ${sim#*.} == 3 ];
then
    cat ${sim}3.mask | $FILET/twoPopnStats_forML $p1 $p2 | python $FILET/normalizeTwoPopnStats.py ${mask}-unmaskfrac $window > ${sim}3.fvout2
else
    cat ${sim}3.mask | $FILET/twoPopnStats_forML $p1 $p2 -c | python $FILET/normalizeTwoPopnStats.py ${mask}-unmaskfrac $window > ${sim}3.fvout2
fi
grep -v "inf" ${sim}3.fvout2 | grep -v "nan" > ${sim}3.fvout
sed -i '/^$/d' ${sim}3.fvout
rm -f ${sim}3.mask
rm -f ${sim}3
rm -f ${sim}3.fvout2