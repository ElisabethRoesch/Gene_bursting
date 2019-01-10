#!/bin/sh

cd /project/home17/er4517/ABC_bursting/ABC_bursting/Code

count=0



n_processes=$1
ABC_file=$2
start_index=$3
end_index=$4
output_file=$5
output_log_file=$6
echo $1
echo "Run started at " `date` > $output_log_file

# module load julia/0.6.0

julia -p $n_processes $ABC_file $start_index $end_index $output_file 2>&1 >> $output_log_file

echo "Run finished at " `date` >> $output_log_file
