
qsub -keo -q long -lnodes=1:ppn=16:cuda10 -v MYARGUMENTS="16 ../../src-test/ABC_SMC_RNAseq_fixedparams.jl 1 16 $ABCDATA/ABC_output/RNAseq/ABCout_RNAseq_lowSF_transformed_keygenes_2i $ABCDATA/ABC_output/log_files/ABC_RNAseq_lowSF_transformed_cuda_2i.log" run_ABC.sh

qsub -keo -q long -lnodes=1:ppn=16:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_1_cell_cycle_2i.jl 1 6 $ABCDATA/ABC_output/RNAseq/lisi_high_var_2i_cellcycle $ABCDATA/ABC_output/log_files/lisi_high_var_2i_cellcycle.log" run_ABC.sh

qsub -keo -q long -lnodes=1:ppn=16:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_1_cell_cycle_2i.jl 1 6 $ABCDATA/ABC_output/RNAseq/lisi_high_var_2i_cellcycle $ABCDATA/ABC_output/log_files/lisi_high_var_2i_cellcycle.log" run_ABC.sh

qsub -keo -q long -lnodes=1:ppn=16:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_1_cell_cycle_2i.jl 1 6 $ABCDATA/ABC_output/RNAseq/lisi_high_var_2i_cellcycle $ABCDATA/ABC_output/log_files/lisi_high_var_2i_cellcycle.log" run_ABC.sh



qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_8_2i2.jl 1 6 ../Data/ABC_output/RNAseq/model_8_2i $ABCDATA/ABC_output/log_files/model_8_2i2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_9_2i2.jl 1 6 ../Data/ABC_output/RNAseq/model_9_2i $ABCDATA/ABC_output/log_files/model_9_2i2.log" run_ABC.sh


#started 10.7.18
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_1_cell_cycle_2i2.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_1_cellcycle_2i_easy $ABCDATA/ABC_output/log_files/lisi_high_var_2i_cellcycle_easy.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_8_2i2.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_8_2i_easy $ABCDATA/ABC_output/log_files/model_8_2i_easy.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_9_2i2.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_2i_easy $ABCDATA/ABC_output/log_files/model_9_2i_easy.log" run_ABC.sh


#started 11.7.18
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_8_2i2.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_8_2i_medium $ABCDATA/ABC_output/log_files/model_8_2i_medium.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_8_2i3.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_8_2i_easy_many $ABCDATA/ABC_output/log_files/model_8_2i_easy_many.log" run_ABC.sh

#started 13.7.18
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_9_2i2.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_2i_medium $ABCDATA/ABC_output/log_files/model_9_2i_medium.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_1_cell_cycle_2i2.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_1_cellcycle_2i_medium $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_medium.log" run_ABC.sh


#started 17.7.18: nanog gene basemodel, basemodel with cellcycle,8,9. Also basemodel 2i medium to compare again.
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_1_serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/model_1_serum_medium $ABCDATA/ABC_output/log_files/model_1_serum_medium.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_1_cell_cycle_serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/model_1_cellcycle_serum_medium $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_medium.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/ABC_SMC_fixedparams_model_9_serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/model_9_serum_medium $ABCDATA/ABC_output/log_files/model_9_serum_medium.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/ABC_SMC_fixedparams_model_8_serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/model_8_serum_medium $ABCDATA/ABC_output/log_files/model_8_serum_medium.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=2:cuda10 -v MYARGUMENTS="2 ../../src-test/ABC_SMC_fixedparams_model_1_2i.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_1_2i_medium $ABCDATA/ABC_output/log_files/model_1_2i_medium.log" run_ABC.sh

#started 17.7.18: nanog gene basemodel, basemodel with cellcycle,8,9. with 2i data.
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/ABC_SMC_fixedparams_model_1_cell_cycle_serum_2.jl 2 2 $ABCDATA/ABC_output/RNAseq/model_1_cellcycle_serum_medium_2 $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_medium_2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/ABC_SMC_fixedparams_model_1_serum_2.jl 2 2 $ABCDATA/ABC_output/RNAseq/model_1_serum_medium_2 $ABCDATA/ABC_output/log_files/model_1_serum_medium_2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/ABC_SMC_fixedparams_model_9_serum_2.jl 2 2 $ABCDATA/ABC_output/RNAseq/model_9_serum_medium_2 $ABCDATA/ABC_output/log_files/model_9_serum_medium_2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/ABC_SMC_fixedparams_model_8_serum_2.jl 2 2 $ABCDATA/ABC_output/RNAseq/model_8_serum_medium_2 $ABCDATA/ABC_output/log_files/model_8_serum_medium_2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/validation_simulation/VALIDATE_easy_model_1.jl 1 1 $ABCDATA/ABC_output/RNAseq/validation_simulation/model_1 $ABCDATA/ABC_output/log_files/VALIDATE_easy_model_1.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/validation_simulation/VALIDATE_easy_model_1.jl 1 1 $ABCDATA/ABC_output/RNAseq/validation_simulation/model_1 $ABCDATA/ABC_output/log_files/VALIDATE_easy_model_1.log" run_ABC.sh
#TODO run high var serum rest
#TODO run high var serum genes in 2i
#TODO run high var 2i genes in serum
#3 are still not finished on 23.7.18: id problem:output folders did not exist !


#started 23.7.18: base model new prior other epsilon
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/testing-priors/model_1_tmax_usual.jl 1 1 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_usual.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/testing-priors/model_1_tmax_large.jl 1 1 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_large.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/ABC_SMC_fixedparams_model_1_serum_2.jl 2 2 $ABCDATA/ABC_output/RNAseq/model_1_serum_medium_2 $ABCDATA/ABC_output/log_files/model_1_serum_medium_2.log" run_ABC.sh

#started 24.7.18
#new runs: baselmodel, epsilons 0.4,0.2,0.1 and n_samples*500 = tmax
#priors: -6,0 for first three and scaling facto 0,20
#try three options, same gene. n_samples=100,200,300 so tmax  50.000, 100.000, 150.000 respectively.
#ids : 87517,87520,87521, prior40:87522,87523,87524. prior 80: 87528/29/30 epsilon min 005: 87531/2/3
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-testing-priors/model_1_tmax_50000.jl 1 1 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_50000.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-testing-priors/model_1_tmax_100000.jl 1 1 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_100000.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-testing-priors/model_1_tmax_150000.jl 1 1 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_150000.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-testing-priors/model_1_tmax_300000.jl 1 1 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_300000_gene1.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-testing-priors/model_1_tmax_300000.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_300000_gene2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-testing-priors/model_1_tmax_300000.jl 3 3 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_300000_gene3.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-testing-priors/model_1_tmax_300000_prior80.jl 1 1 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_300000_gene1_prior80.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-testing-priors/model_1_tmax_300000_prior80.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_300000_gene2_prior80.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-testing-priors/model_1_tmax_300000_prior80.jl 3 3 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_300000_gene3_prior80.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-testing-priors/model_1_tmax_300000_005.jl 1 1 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_300000_gene1_005.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-testing-priors/model_1_tmax_300000_005.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_300000_gene2_005.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-testing-priors/model_1_tmax_300000_005.jl 3 3 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_tmax_300000_gene3_005.log" run_ABC.sh



#last thee from yesterday are still running
#starting 25.7: ?,87539
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors/model-1-100-fix.jl 1 1 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_100.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors/model-1-200-fix.jl 1 1 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_200.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors/model-1-500-fix.jl 1 1 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_500.log" run_ABC.sh


#26.7:87548,87546,87547
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors/model-1-100-fix.jl 1 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_100_50000.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors/model-1-200-fix.jl 1 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_200_100000.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors/model-1-500-fix.jl 1 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1 $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_500_250000.log" run_ABC.sh


#27.7
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors/model-1-100-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_100_50000_2i.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors/model-1-200-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_200_100000_2i.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors/model-1-500-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_500_250000_2i.log" run_ABC.sh

qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors/model-1-100-fix-serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/serum $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_100_50000_serum.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors/model-1-200-fix-serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/serum $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_200_100000_serum.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors/model-1-500-fix-serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/serum $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_500_250000_serum.log" run_ABC.sh


qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-1000/model-1-100-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/1000/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_100_50000_2i_1000.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-1000/model-1-200-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/1000/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_200_100000_2i_1000.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-1000/model-1-500-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/1000/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_500_250000_2i_1000.log" run_ABC.sh

qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-1000/model-1-100-fix-serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/1000/serum $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_100_50000_serum_1000.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-1000/model-1-200-fix-serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/1000/serum $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_200_100000_serum_1000.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-1000/model-1-500-fix-serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/1000/serum $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_500_250000_serum_1000.log" run_ABC.sh

#30.7
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-500/model-1-100-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/nanog/500/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_100_50000_2i.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-500/model-1-200-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/nanog/500/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_200_100000_2i.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-1000/model-1-100-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/nanog/1000/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_100_50000_2i_1000.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-1000/model-1-200-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/nanog/1000/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_200_100000_2i_1000.log" run_ABC.sh


qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-250/model-1-100-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/nanog/200/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_100_25000_2i.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-250/model-1-200-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/nanog/200/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_200_50000_2i.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-250/model-1-500-fix.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/nanog/200/2i $ABCDATA/ABC_output/log_files/model_1_2i_0.1_fix_500_125000_2i_1000.log" run_ABC.sh

qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-250/model-1-100-fix-serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/nanog/200/serum $ABCDATA/ABC_output/log_files/model_1_serum_0.1_fix_100_25000_2i_200.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-250/model-1-200-fix-serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/nanog/200/serum $ABCDATA/ABC_output/log_files/model_1_serum_0.1_fix_200_50000_2i_200.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-250/model-1-500-fix-serum.jl 2 2 $ABCDATA/ABC_output/RNAseq/testing_priors/model_1/nanog/200/serum $ABCDATA/ABC_output/log_files/model_1_serum_0.1_fix_500_125000_serum_200.log" run_ABC.sh




#1.8
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-serum.jl 6 6 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/serum_bimodule $ABCDATA/ABC_output/log_files/model_1_serum_bimodule_500_250000_2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-serum.jl 8 8 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/serum_bimodule $ABCDATA/ABC_output/log_files/model_1_serum_bimodule_500_250000_1.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-serum.jl 10 10 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/serum_bimodule $ABCDATA/ABC_output/log_files/model_1_serum_bimodule_500_250000_3.log" run_ABC.sh

qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-2i.jl 4 4 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/2i_rightpeak $ABCDATA/ABC_output/log_files/model_1_2i_rightpeak_500_250000_2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-2i.jl 10 10 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/2i_rightpeak $ABCDATA/ABC_output/log_files/model_1_2i_rightpeak_500_250000_3.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-2i.jl 11 11 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/2i_rightpeak $ABCDATA/ABC_output/log_files/model_1_2i_rightpeak_500_250000_1.log" run_ABC.sh


#2.8
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-2i-genes-2i-data.jl 1 6 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/all_high_var_in_2i_genes_in_2i $ABCDATA/ABC_output/log_files/model_1_all_high_var_in_2i_genes_in_2i.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-2i-genes-serum-data.jl 1 6 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/all_high_var_in_2i_genes_in_serum $ABCDATA/ABC_output/log_files/model_1_all_high_var_in_2i_genes_in_serum.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=11:node030 -v MYARGUMENTS="11 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-serum-genes-2i-data.jl 1 11 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/all_high_var_in_serum_genes_in_2i $ABCDATA/ABC_output/log_files/model_1_all_high_var_in_serum_genes_in_2i.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=11:node030 -v MYARGUMENTS="11 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-serum-genes-serum-data.jl 1 11 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/all_high_var_in_serum_genes_in_serum $ABCDATA/ABC_output/log_files/model_1_all_high_var_in_serum_genes_in_serum.log" run_ABC.sh


#3.8 cell cycle nanog t1=t2 = 12500/25000
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_serum_genes_2i_data_12500_12500.jl 2 2 $ABCDATA/ABC_output/RNAseq/cell_cycle/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_genes_2i_data_12500_12500.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_serum_genes_2i_data_25000_25000.jl 2 2 $ABCDATA/ABC_output/RNAseq/cell_cycle/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_genes_2i_data_25000_25000.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_serum_genes_serum_data_12500_12500.jl 2 2 $ABCDATA/ABC_output/RNAseq/cell_cycle/serum_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_genes_serum_data_12500_12500.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_serum_genes_serum_data_25000_25000.jl 2 2 $ABCDATA/ABC_output/RNAseq/cell_cycle/serum_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_genes_serum_data_25000_25000.log" run_ABC.sh



#6.8 cell cycle nanog t1!=t2  start from home!
#CCna2
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_2i_data_8333_16666.jl 3 3 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_2i_data_8333_16666_CCna2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_2i_data_16666_33332.jl 3 3 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_2i_data_16666_33332_CCna2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_serum_data_8333_16666.jl 3 3 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_serum_data_8333_16666_CCna2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_serum_data_16666_33332.jl 3 3 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_serum_data_16666_33332_CCna2.log" run_ABC.sh
#CCnb1
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_2i_data_8333_16666.jl 5 5 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_2i_data_8333_16666_CCnb1.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_2i_data_16666_33332.jl 5 5 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_2i_data_16666_33332_CCnb1.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_serum_data_8333_16666.jl 5 5 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_serum_data_8333_16666_CCnb1.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_serum_data_16666_33332.jl 5 5 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_serum_data_16666_33332_CCnb1.log" run_ABC.sh
#Hes1
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_serum_genes_2i_data_8333_16666.jl 7 7 $ABCDATA/ABC_output/RNAseq/cell_cycle/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_genes_2i_data_8333_16666_Hes1.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_serum_genes_2i_data_16666_33332.jl 7 7 $ABCDATA/ABC_output/RNAseq/cell_cycle/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_genes_2i_data_16666_33332_Hes1.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_serum_genes_serum_data_8333_16666.jl 7 7 $ABCDATA/ABC_output/RNAseq/cell_cycle/serum_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_genes_serum_data_8333_16666_Hes1.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_serum_genes_serum_data_16666_33332.jl 7 7 $ABCDATA/ABC_output/RNAseq/cell_cycle/serum_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_genes_serum_data_16666_33332_Hes1.log" run_ABC.sh
#basemodel all genes
qsub -keo -q long -lnodes=1:ppn=11:cpu:node030 -v MYARGUMENTS="11 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-serum-genes-2i-data.jl 1 11 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/all_high_var_in_serum_genes_in_2i $ABCDATA/ABC_output/log_files/model_1_all_high_var_in_serum_genes_in_2i.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=11:cpu:node030 -v MYARGUMENTS="11 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-serum-genes-serum-data.jl 1 11 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/all_high_var_in_serum_genes_in_serum $ABCDATA/ABC_output/log_files/model_1_all_high_var_in_serum_genes_in_serum.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=4:cpu:node030 -v MYARGUMENTS="4 ../../src-test/abc-smc-fix-priors-500/model-1-500-fix-2i-genes-serum-data.jl 3 6 $ABCDATA/ABC_output/RNAseq/explore_distr/model_1/all_high_var_in_2i_genes_in_serum $ABCDATA/ABC_output/log_files/model_1_all_high_var_in_2i_genes_in_serum.log" run_ABC.sh



#11.8
#CCna2
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_2i_data_8333_16666.jl 2 2 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_2i_data_8333_16666_CCna2_real.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_2i_data_16666_33332.jl 2 2 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_2i_data_16666_33332_CCna2_real.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_serum_data_8333_16666.jl 2 2 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_serum_data_8333_16666_CCna2_real.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_serum_data_16666_33332.jl 2 2 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_serum_data_16666_33332_CCna2_real.log" run_ABC.sh
#CCnb1
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_2i_data_8333_16666.jl 6 6 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_2i_data_8333_16666_CCnb1_real.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_2i_data_16666_33332.jl 6 6 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_2i_data_16666_33332_CCnb1_real.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_serum_data_8333_16666.jl 6 6 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_serum_data_8333_16666_CCnb1_real.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_serum_data_16666_33332.jl 6 6 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_serum_data_16666_33332_CCnb1_real.log" run_ABC.sh
#test simulated data
qsub -keo -q long -lnodes=1:ppn=1:cuda10 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate-base-model.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_1 $ABCDATA/ABC_output/log_files/model_1_validate_295.log" run_ABC.sh


#all high var auf 5 times cc
qsub -keo -q long -lnodes=1:ppn=6:cpu:node030 -v MYARGUMENTS="6 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_2i_data_8333_16666.jl 1 6 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_2i_data_8333_16666_all_real.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:cpu:node030 -v MYARGUMENTS="6 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_serum_data_8333_16666.jl 1 6 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_serum_data_8333_16666_all_real.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=11:cpu:node030 -v MYARGUMENTS="11 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_serum_genes_2i_data_8333_16666.jl 1 11 $ABCDATA/ABC_output/RNAseq/cell_cycle/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_genes_2i_data_8333_16666_all.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=11:cpu:node030 -v MYARGUMENTS="11 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_serum_genes_serum_data_8333_16666.jl 1 11 $ABCDATA/ABC_output/RNAseq/cell_cycle/serum_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_genes_serum_data_8333_16666_all.log" run_ABC.sh


#todo########################### this is same as above but 10 swiches
qsub -keo -q long -lnodes=1:ppn=6:cpu:node030 -v MYARGUMENTS="6 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_2i_data_16666_33332.jl 1 6 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_2i_data_16666_33332_all_real.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:cpu:node030 -v MYARGUMENTS="6 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_2i_genes_serum_data_16666_33332.jl 1 6 $ABCDATA/ABC_output/RNAseq/cell_cycle/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_2i_genes_serum_data_16666_33332_all_real.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=11:cpu:node030 -v MYARGUMENTS="11 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_serum_genes_2i_data_16666_33332.jl 1 11 $ABCDATA/ABC_output/RNAseq/cell_cycle/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_genes_2i_data_16666_33332_all.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=11:cpu:node031 -v MYARGUMENTS="11 ../../src-test/cell_cycle/ABC_SMC_fixedparams_model_1_cell_cycle_serum_genes_serum_data_16666_33332.jl 1 11 $ABCDATA/ABC_output/RNAseq/cell_cycle/serum_genes/serum_data $ABCDATA/ABC_output/log_files/model_1_cellcycle_serum_genes_serum_data_16666_33332_all.log" run_ABC.sh
##############################

#12.8
qsub -keo -q long -lnodes=1:ppn=5:cpu:node031 -v MYARGUMENTS="5 ../../src-test/model_8/ABC_SMC_model_8_2i_genes_2i_data.jl 1 5 $ABCDATA/ABC_output/RNAseq/model_8/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_8_2i_genes_2i_data.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=5:cpu:node031 -v MYARGUMENTS="5 ../../src-test/model_8/ABC_SMC_model_8_2i_genes_serum_data.jl 1 5 $ABCDATA/ABC_output/RNAseq/model_8/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_8_2i_genes_serum_data.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=5:cpu:node031 -v MYARGUMENTS="5 ../../src-test/model_8/ABC_SMC_model_8_serum_genes_serum_data.jl 1 5 $ABCDATA/ABC_output/RNAseq/model_8/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_8_serum_genes_serum_data.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=5:cpu:node031 -v MYARGUMENTS="5 ../../src-test/model_8/ABC_SMC_model_8_serum_genes_2i_data.jl 1 5 $ABCDATA/ABC_output/RNAseq/model_8/serum_genes/serum_data $ABCDATA/ABC_output/log_files/model_8_serum_genes_2i_data.log" run_ABC.sh


#sam but 0.4 and broader prior
qsub -keo -q long -lnodes=1:ppn=5:cpu:node031 -v MYARGUMENTS="5 ../../src-test/model_8/ABC_SMC_model_8_2i_genes_2i_data.jl 1 5 $ABCDATA/ABC_output/RNAseq/model_8_0.4_5/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_8_2i_genes_2i_data3.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=5:cpu:node031 -v MYARGUMENTS="5 ../../src-test/model_8/ABC_SMC_model_8_2i_genes_serum_data.jl 1 5 $ABCDATA/ABC_output/RNAseq/model_8_0.4_5/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_8_2i_genes_serum_data3.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=5:cpu:node031 -v MYARGUMENTS="5 ../../src-test/model_8/ABC_SMC_model_8_serum_genes_serum_data.jl 1 5 $ABCDATA/ABC_output/RNAseq/model_8_0.4_5/serum_genes/serum_data  $ABCDATA/ABC_output/log_files/model_8_serum_genes_serum_data3.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=5:cpu:node031 -v MYARGUMENTS="5 ../../src-test/model_8/ABC_SMC_model_8_serum_genes_2i_data.jl 1 5 $ABCDATA/ABC_output/RNAseq/model_8_0.4_5/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_8_serum_genes_2i_data3.log" run_ABC.sh

#9model
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_2i_genes_2i_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_0.4_5_5/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_9_2i_genes_2i_data.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:cuda10 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_2i_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_0.4_5_5/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_9_2i_genes_serum_data.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=10:cuda10 -v MYARGUMENTS="10 ../../src-test/model_9/ABC_SMC_model_9_serum_genes_serum_data.jl 1 10 $ABCDATA/ABC_output/RNAseq/model_9_0.4_5_5/serum_genes/serum_data  $ABCDATA/ABC_output/log_files/model_9_serum_genes_serum_data.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=10:cuda10 -v MYARGUMENTS="10 ../../src-test/model_9/ABC_SMC_model_9_serum_genes_2i_data.jl 1 10 $ABCDATA/ABC_output/RNAseq/model_9_0.4_5_5/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_9_serum_genes_2i_data.log" run_ABC.sh



# 14.8
#87731
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate-base-model.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_1 $ABCDATA/ABC_output/log_files/model_1_validate_295.log" run_ABC.sh
#87732-87735
qsub -keo -q long -lnodes=1:ppn=6:node030 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_2i_genes_2i_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_20/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_9_2i_genes_2i_data_20.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node030 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_2i_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_20/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_9_2i_genes_serum_data_20.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node030 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_serum_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_20/serum_genes/serum_data  $ABCDATA/ABC_output/log_files/model_9_serum_genes_serum_data_20.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node030 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_serum_genes_2i_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_20/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_9_serum_genes_2i_data_20.log" run_ABC.sh
#87737.
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate-base-model0.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_1 $ABCDATA/ABC_output/log_files/model_1_validate_295.log" run_ABC.sh



#15.8: mini epsilon 0.2
qsub -keo -q long -lnodes=1:ppn=6:node030 -v MYARGUMENTS="6 ../../src-test/model_9_02/ABC_SMC_model_9_2i_genes_2i_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_david_02/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_9_2i_genes_2i_data_david2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node030 -v MYARGUMENTS="6 ../../src-test/model_9_02/ABC_SMC_model_9_2i_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_david_02/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_9_2i_genes_serum_data_david2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node030 -v MYARGUMENTS="6 ../../src-test/model_9_02/ABC_SMC_model_9_serum_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_david_02/serum_genes/serum_data  $ABCDATA/ABC_output/log_files/model_9_serum_genes_serum_data_david2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node030 -v MYARGUMENTS="6 ../../src-test/model_9_02/ABC_SMC_model_9_serum_genes_2i_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_david_02/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_9_serum_genes_2i_data_david2.log" run_ABC.sh
#mini epsilon 0.4
qsub -keo -q long -lnodes=1:ppn=6:node030 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_2i_genes_2i_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_david/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_9_2i_genes_2i_data_david.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node031 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_2i_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_david/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_9_2i_genes_serum_data_david.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node031 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_serum_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_david/serum_genes/serum_data  $ABCDATA/ABC_output/log_files/model_9_serum_genes_serum_data_david.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node031 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_serum_genes_2i_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_david/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_9_serum_genes_2i_data_david.log" run_ABC.sh
#



#simulated data for m8,m9
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_light_fb.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/model_9_validate_light_fb.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_heavy_fb.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/model_9_validate_heavy_fb.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_8_light_fb.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_8 $ABCDATA/ABC_output/log_files/model_8_validate_light_fb.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_8_heavy_fb.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_8 $ABCDATA/ABC_output/log_files/model_8_validate_heavy_fb.log" run_ABC.sh

#model 8:  mini epsilon 0.2
qsub -keo -q long -lnodes=1:ppn=5:cpu:node031 -v MYARGUMENTS="6 ../../src-test/model_8/ABC_SMC_model_8_2i_genes_2i_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_8_02/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_8_2i_genes_2i_data3_02.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=5:cpu:node031 -v MYARGUMENTS="6 ../../src-test/model_8/ABC_SMC_model_8_2i_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_8_02/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_8_2i_genes_serum_data3_02.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=5:cpu:node031 -v MYARGUMENTS="6 ../../src-test/model_8/ABC_SMC_model_8_serum_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_8_02/serum_genes/serum_data  $ABCDATA/ABC_output/log_files/model_8_serum_genes_serum_data3_02.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=5:cpu:node031 -v MYARGUMENTS="6 ../../src-test/model_8/ABC_SMC_model_8_serum_genes_2i_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_8_02/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_8_serum_genes_2i_data3.log_02" run_ABC.sh

############## 21.8. nach david skype ##################
#model 9 new priors.
qsub -keo -q long -lnodes=1:ppn=4:node030 -v MYARGUMENTS="4 ../../src-test/model_9_02/ABC_SMC_model_9_2i_genes_2i_data.jl 1 4 $ABCDATA/ABC_output/RNAseq/model_9_huge/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_9_2i_genes_2i_data_huge.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=4:node030 -v MYARGUMENTS="4 ../../src-test/model_9_02/ABC_SMC_model_9_2i_genes_serum_data.jl 1 4 $ABCDATA/ABC_output/RNAseq/model_9_huge/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_9_2i_genes_serum_data_huge.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node030 -v MYARGUMENTS="6 ../../src-test/model_9_02/ABC_SMC_model_9_serum_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_huge/serum_genes/serum_data  $ABCDATA/ABC_output/log_files/model_9_serum_genes_serum_data_huge.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node030 -v MYARGUMENTS="6 ../../src-test/model_9_02/ABC_SMC_model_9_serum_genes_2i_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_huge/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_9_serum_genes_2i_data_huge.log" run_ABC.sh

#model 9 new priors no deact.
qsub -keo -q long -lnodes=1:ppn=4:node031 -v MYARGUMENTS="4 ../../src-test/model_9_02_no_deact/ABC_SMC_model_9_2i_genes_2i_data.jl 1 4 $ABCDATA/ABC_output/RNAseq/model_9_huge_no_deact/2i_genes/2i_data $ABCDATA/ABC_output/log_files/model_9_2i_genes_2i_data_huge_no_deact.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=4:node031 -v MYARGUMENTS="4 ../../src-test/model_9_02_no_deact/ABC_SMC_model_9_2i_genes_serum_data.jl 1 4 $ABCDATA/ABC_output/RNAseq/model_9_huge_no_deact/2i_genes/serum_data $ABCDATA/ABC_output/log_files/model_9_2i_genes_serum_data_huge_no_deact.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node031 -v MYARGUMENTS="6 ../../src-test/model_9_02_no_deact/ABC_SMC_model_9_serum_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_huge_no_deact/serum_genes/serum_data  $ABCDATA/ABC_output/log_files/model_9_serum_genes_serum_data_huge_no_deact.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node031 -v MYARGUMENTS="6 ../../src-test/model_9_02_no_deact/ABC_SMC_model_9_serum_genes_2i_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_huge_no_deact/serum_genes/2i_data $ABCDATA/ABC_output/log_files/model_9_serum_genes_2i_data_huge_no_deact.log" run_ABC.sh

#simulated data model 9 new priors with and without deact
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_light_fb_01_no_deact.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_light_fb_01_no_deact.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_light_fb_02_no_deact.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_light_fb_02_no_deact.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_light_fb_01.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_light_fb_01.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_light_fb_02.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_light_fb_02.log" run_ABC.sh



######tursday late at night start davids suggestions83,84
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_100_02.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_02.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_100_01.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_01.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_500_01.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_500_01.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_500_02.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_500_02.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_100_01_250.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_01_250.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_100_02_250.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_02_250.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_500_01_250.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_500_01_250.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_500_02_250.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_500_02_250.log" run_ABC.sh

#friday: model 9 with more precise singe cell data. we change from 250 to 500 87791,2,3,4
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_100_01_500.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_01_500.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_100_02_500.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_02_500.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_500_01_500.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_500_01_500.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node031 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_500_02_500.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_500_02_500.log" run_ABC.sh

#friday: shame on you: k was 10^2 instead on 1063 as david said. 87795,6,7,8
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_100_01_500_3.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_01_500_3.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_100_02_500_3.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_02_500_3.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_500_01_500_3.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_500_01_500_3.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_500_02_500_3.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_500_02_500_3.log" run_ABC.sh

#friday: model 8 david kind of obs: 87799,00,01,02
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_8_100_01_250.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_8_100_01_250.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_8_100_02_250.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_8_100_02_250.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_8_100_01_500.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_8_100_01_500.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_8_100_02_500.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_8_100_02_500.log" run_ABC.sh


#friday: cc base model kind  87803,04,05,06
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_cc_100_01_295.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_cc_100_01_295.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_cc_100_02_295.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_cc_100_02_295.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_cc_500_01_295.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_cc_500_01_295.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_cc_500_02_295.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_cc_500_02_295.log" run_ABC.sh

#friday: cc only one as ergodic 87807,8
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_cc_100_01_295_only_one.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_cc_100_01_295_only_one.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_cc_100_02_295_only_one.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_cc_100_02_295_only_one.log" run_ABC.sh


#saturday: model 9: newTrue=log10.([0.001, 0, 0.0001, 0.01, 1000.])
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_100_02_500_0c2.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_02_500_0c2.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/validate_model_9_david_combi_100_02_500_3_0c2.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9 $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_02_500_3_0c2.log" run_ABC.sh



qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/no-deact-mo/validate_model_9_david_combi_100_02_250_3_0c2_new_prior.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9/new $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_02_250_3_0c2_new_prior.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/no-deact-mo/validate_model_9_david_combi_100_02_250_3_0c2_new_prior_save.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9/new $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_02_250_3_0c2_new_prior_save.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/no-deact-mo/validate_model_9_david_combi_100_01_250_3_0c2_new_prior_save.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9/new $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_01_250_3_0c2_new_prior_save.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/no-deact-k/validate_model_9_david_combi_100_01_250_3_0c2_no_k.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9/new $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_01_250_3_0c2_no_k.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/no-deact-k/validate_model_9_david_combi_100_02_250_3_0c2_no_k.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9/new $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_02_250_3_0c2_no_k.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=1:node030 -v MYARGUMENTS="1 ../../src-test/abc-smc-validation-simulation/no-deact-k/validate_model_9_david_combi_100_02_250_3_0c2_no_k_save.jl 1 1  $ABCDATA/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9/new $ABCDATA/ABC_output/log_files/validate_model_9_david_combi_100_02_250_3_0c2_no_k_save.log" run_ABC.sh


validate_model_9_david_combi_100_01_250_3_0c2_new_prior_save
validate_model_9_david_combi_100_01_250_3_0c2_no_k
validate_model_9_david_combi_100_02_250_3_0c2_no_k_save




qsub -keo -q long -lnodes=1:ppn=6:node031 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_serum_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_david/serum_genes/serum_data  $ABCDATA/ABC_output/log_files/model_9_serum_genes_serum_data_david_n.log" run_ABC.sh
qsub -keo -q long -lnodes=1:ppn=6:node031 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_serum_genes_serum_data0.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_david/serum_genes/serum_data  $ABCDATA/ABC_output/log_files/model_9_serum_genes_serum_data_david_n0.log" run_ABC.sh


qsub -keo -q long -lnodes=1:ppn=6:node031 -v MYARGUMENTS="6 ../../src-test/model_9/ABC_SMC_model_9_serum_genes_serum_data.jl 1 6 $ABCDATA/ABC_output/RNAseq/model_9_david/serum_genes/serum_data  $ABCDATA/ABC_output/log_files/model_9_serum_genes_serum_data_david_n.log" run_ABC.sh
