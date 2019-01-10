using PyCall, PyPlot; @pyimport seaborn as sns
using ABC_Bursting

names=["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]

for name in names
  pathtofile_lisi = string(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_1_cellcycle_2i_easy/model_1_cellcycle_easy_",name),"_2i_0.8_100.txt")
  pathtoplot = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/model-1-cellcycle-easy/",name)
  plot_one_abacus_file_m1(pathtofile_lisi, pathtoplot, 3, 100, 4)
end
