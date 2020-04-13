#!/bin/bash

while read each_subject
do

# copy the batch script template
cp /work/TEDDY/script/summary/error_correction.batch  $each_subject/${each_subject}_correct.batch

# change 'SUBJECTID' in the template into your subject id
sed -i "s/SUBJECTID/$each_subject/g"  $each_subject/${each_subject}_correct.batch

# set the parameters of batch script
#sed -i "s/omicsbio/normal/g" $each_subject/${each_subject}_correct.batch
#sed  -i "s/time=72/time=48/g"   $each_subject/${each_subject}_correct.batch

# submit batch file into server
sbatch  $each_subject/${each_subject}_correct.batch

done < subjectID_list.txt