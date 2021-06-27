#!/bin/bash

while read i 
do

cp /work/TEDDY/script/summary/abundance_binning/merge_abunance/merge.batch /work/TEDDY/abundance/finished/$i/${i}_merge.batch
sed -i "s/SUBJECTID/$i/g" /work/TEDDY/abundance/finished/$i/${i}_merge.batch
sbatch /work/TEDDY/abundance/finished/$i/${i}_merge.batch


done < subjectID_list.txt
