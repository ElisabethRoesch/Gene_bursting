using PyCall, PyPlot; @pyimport seaborn as sns
using PyPlot
using ABC_Bursting
# using PyCall, PyPlot; @pyimport matplotlib.pyplot as plat

genes_2i_high_var = ["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]
genes_serum_high_var = ["ENSMUSG00000003032", "ENSMUSG00000012396", "ENSMUSG00000018604", "ENSMUSG00000020717", "ENSMUSG00000021255", "ENSMUSG00000021835", "ENSMUSG00000022528", "ENSMUSG00000025056", "ENSMUSG00000038793", "ENSMUSG00000048402", "ENSMUSG00000051176"]
data_2i = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_2i_filtereddata.txt")
data_serum = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_serum_filtereddata.txt")
binsize =20
for i in 1:length(genes_2i_high_var)
    gene = genes_serum_high_var[i]
    #read and store experimental data_2i for gene "gene"
    one_gene_experimental_result_2i = Get_single_gene_RNAseq_data(gene,data_2i)
    one_gene_experimental_result_serum = Get_single_gene_RNAseq_data(gene,data_serum)
    #plot experimental data
    fig,ax = PyPlot.subplots()
    sns.distplot(one_gene_experimental_result_2i,norm_hist=true, bins=binsize, kde= true, hist=true,rug= false ,  color = "red", label = "experimental 2i",axlabel = "Product abundance",ax=ax)
    sns.distplot(one_gene_experimental_result_serum,norm_hist=true, bins=binsize, kde= true, hist=true,rug= false ,  color = "blue", label = "experimental serum",axlabel = "Product abundance",ax=ax)
    name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/real_data/genes_high_var_2i/", gene,".png")
    savefig(name)
end
