
#MYARGUMENTS: Modify as required
# arg 1 = n_processes
# arg 2 = ABC_file
# arg 3 = start_index
# arg 4 = end_index
# arg 5 = output_folder
# arg 6 = output_log_file

#locate to code file
cd /project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Code/src/sh
#make the output folder
mkdir ./output_files/ABC_RNAseq_sampled

#Run ABC on the servers:
#Will run genes 1:16 on cuda10 using 16 processes
qsub -keo -q long -lnodes=1:ppn=16:cuda10 -v MYARGUMENTS="16 ./ABC_SMC_RNAseq_fixedparams.jl 1 16 ../Data/ABC_output/RNAseq/ABCout_RNAseq_lowSF_transformed_keygenes_2i ../Data/ABC_output/log_files/ABC_RNAseq_lowSF_transformed_cuda_2i.log" run_ABC.sh

#Will run genes 17:43 on node030 using 24 processes
qsub -keo -q long -lnodes=1:ppn=24:node030 -v MYARGUMENTS="24 ./ABC_SMC_RNAseq_fixedparams.jl 17 43 ../Data/ABC_output/RNAseq/ABCout_RNAseq_lowSF_transformed_keygenes_2i ../Data/ABC_output/log_files/ABC_RNAseq_lowSF_transformed_node030_2i.log" run_ABC.sh

#Will run genes 44:68 on node031 using 24 processes
qsub -keo -q long -lnodes=1:ppn=24:node031 -v MYARGUMENTS="24 ./ABC_SMC_RNAseq_fixedparams.jl 44 68 ../Data/ABC_output/RNAseq/ABCout_RNAseq_lowSF_transformed_keygenes_2i ../Data/ABC_output/log_files/ABC_RNAseq_lowSF_transformed_node031_2i.log" run_ABC.sh


#Outputs will be saved in the folder ../Data/ABC_output/RNAseq/ABCout_RNAseq_lowSF_transformed_keygenes_2i
#Logs will be saved in ../Data/ABC_output/log_files/ABC_RNAseq_lowSF_transformed_node031_2i.log
