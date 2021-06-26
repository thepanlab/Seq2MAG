#!/bin/bash

cd /work/TEDDY/cluster/kofam/script

while read i
do
cp /work/TEDDY/script/summary/kofam/kofam.batch ${i}_kofam.batch

sed -i "s/SUBJECTID/$i/g"  ${i}_kofam.batch
#sed -i "s/normal/omicsbio/g" ${i}_kofam.batch
#sed  -i "s/time=48/time=72/g"   ${i}_kofam.batch

sbatch  ${i}_kofam.batch

done < subjectID_list.txt

