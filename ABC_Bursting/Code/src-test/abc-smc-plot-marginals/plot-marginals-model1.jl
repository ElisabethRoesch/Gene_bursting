using PyCall, PyPlot; @pyimport seaborn as sns
using ABC_Bursting

names=["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]

names=["ENSMUSG00000006398","ENSMUSG00000027715", "ENSMUSG00000029177"]
names=["ENSMUSG00000006398"]

for name in names
  #_300.000_2i_0.1_600
  fileending="_300.000_2i_0.1_prior80_600"
  # pathtofile_lisi = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/model1_",name,"_40.000_2i_0.1_100.txt")
  pathtofile_lisi = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/model1_",name,fileending,".txt")
  pathtoplot = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/testing-priors/2_6_40_80/marginals/",name,fileending,"_marginal.png")
  plot_one_abacus_file_m1(pathtofile_lisi, pathtoplot, 3, 600, 4)
end
