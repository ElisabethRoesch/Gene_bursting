#!/bin/sh

cd /cluster/home/kw1016/Documents/Project_3/Julia_code/Code

#Define the locationg and name of the new log file after >
echo "Run started at " `date` > "/cluster/home/kw1016/Documents/Project_3/Julia_code/Data/ABC_output/log_files/feedback_89_3_nanog.log"

module load julia/0.6.0

#Define the ABC run details
  #number of processes e.g. 4
  #Name of the ABC file to run e.g. "./ABC_SMC_RNAseq_nanog_fixedparams.jl"
  #Name of log file
julia -p 4 "./ABC_SMC_RNAseq_nanog_fixedparams.jl" 2>&1 >> "/cluster/home/kw1016/Documents/Project_3/Julia_code/Data/ABC_output/log_files/feedback_89_3_nanog.log"

echo "Run finished at " `date` >> "/cluster/home/kw1016/Documents/Project_3/Julia_code/Data/ABC_output/log_files/feedback_89_3_nanog.log"


#Run using:
#qsub -keo -q long -lnodes=1:ppn=1:cuda10 run_singleABC.sh
