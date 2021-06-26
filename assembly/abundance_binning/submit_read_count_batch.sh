#!/bin/bash

cd /work/TEDDY/abundance
while read i
do
mkdir $i
cp /work/TEDDY/assembly/metaspades/${i}/scaffolds.fasta  $i/${i}_scaffolds.fasta
cp /work/TEDDY/script/summary/abundance_binning/read_count_new.batch $i/${i}_read_count.batch

sed -i "s/SUBJECTID/$i/g"  $i/${i}_read_count.batch
#sed -i "s/omicsbio/normal/g" $i/${i}_read_count.batch
#sed  -i "s/time=48/time=72/g"   $i/${i}_read_count.batch
sbatch  $i/${i}_read_count.batch

done < /work/TEDDY/script/summary/abundance_binning/subjectID_list.txt

