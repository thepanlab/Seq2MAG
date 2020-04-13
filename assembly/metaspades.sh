#!/bin/bash


while read i
do

# copy the batch script template
cp /work/TEDDY/assembly/metaspades_model.batch  $i/${i}_metaspades.batch

# change 'SUBJECTID' in the template into your subject id
sed -i "s/SUBJECTID/$i/g"  $i/${i}_metaspades.batch

# set the parameters of batch script
sed -i "s/omicsbio/normal/g" $i/${i}_metaspades.batch
#sed  -i "s/time=48/time=72/g"   $i/${i}_metaspades.batch

# submit batch file into server
sbatch  $i/${i}_metaspades.batch

done < subjectID_list.txt
