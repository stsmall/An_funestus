#!bin/bash

#$1 path to reads
#$2 base name

#assembly w2rap-contigger
for k in 200 225 260 275 300;do
/afs/crc.nd.edu/user/s/ssmall2/bin/01_unipaths -o ${1}/w2rap_EC-${k}-July62016 -p ${2} -r ${1}*.1.fastq,${1}*.2.fastq -t 40
/afs/crc.nd.edu/user/s/ssmall2/bin/02_qgraph -o ${1}/w2rap_EC-${k}-July62016 -p ${2} -K $k -t 40
/afs/crc.nd.edu/user/s/ssmall2/bin/03_clean -o ${1}/w2rap_EC-${k}-July62016 -p ${2} -t 40
/afs/crc.nd.edu/user/s/ssmall2/bin/04_patching -o ${1}/w2rap_EC-${k}-July62016 -p ${2} -t 40
/afs/crc.nd.edu/user/s/ssmall2/bin/05_simplify -o ${1}/w2rap_EC-${k}-July62016 -p ${2} -t 40
/afs/crc.nd.edu/user/s/ssmall2/bin/06_scaffolding -o ${1}/w2rap_EC-${k}-July62016 -p ${2} -t 40
done