using ABC_Bursting,PyCall, PyPlot
@pyimport seaborn as sns
@pyimport pandas as pd


gene_names = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/gene_names/new_genes.txt",'\t','\n')
gene_names_en = gene_names[1:end,:,:]
real_names = Dict()
gene_names_en[2,3]
for gene_counter in 1:length(gene_names_en[:,1,1])
   real_names[gene_names_en[gene_counter]] = gene_names_en[gene_counter,2]
end
gene_names = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/gene_names/new_genes.txt",'\t','\n')
gene_names_en = gene_names[1:end,:,:]
gene_names_en[2,3]
for gene_counter in 1:length(gene_names_en[:,1,1])
   real_names[gene_names_en[gene_counter]] = gene_names_en[gene_counter,2]
end
genes_2i_high_var = ["ENSMUSG00000006589", "ENSMUSG00000019539", "ENSMUSG00000020571", "ENSMUSG00000020954", "ENSMUSG00000022881", "ENSMUSG00000024646", "ENSMUSG00000025130", "ENSMUSG00000025359", "ENSMUSG00000047751", "ENSMUSG00000050953", "ENSMUSG00000062270"]

# genes_2i_high_var = ["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]
# genes_serum_high_var = ["ENSMUSG00000003032", "ENSMUSG00000012396", "ENSMUSG00000018604", "ENSMUSG00000020717", "ENSMUSG00000021255", "ENSMUSG00000021835", "ENSMUSG00000022528", "ENSMUSG00000025056", "ENSMUSG00000038793", "ENSMUSG00000048402", "ENSMUSG00000051176"]
data_2i = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_2i_filtereddata.txt")
data_serum = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_serum_filtereddata.txt")
binsize =20

plot_cs = []
plot_conds =[]
plot_names =[]

#here
for i in 1:length(genes_2i_high_var)
   #here
    gene = genes_2i_high_var[i]
    one_gene_experimental_result_2i = Get_single_gene_RNAseq_data(gene,data_2i)
    one_gene_experimental_result_serum = Get_single_gene_RNAseq_data(gene,data_serum)
       n_acc_2i =length(one_gene_experimental_result_2i)
    n_acc_serum =length(one_gene_experimental_result_serum)
    #here
    my_label=string(real_names[genes_2i_high_var[i]])
    c =  vcat(one_gene_experimental_result_2i,one_gene_experimental_result_serum)
    x1 = fill("2i",n_acc_2i)
    x2 = fill("serum",n_acc_serum)
    cond = vcat(x1,x2)
    names = fill(my_label,n_acc_2i+n_acc_serum)
    plot_cs=vcat(plot_cs,c)
    plot_conds=vcat(plot_conds,cond)
    plot_names=vcat(plot_names,names)
end

gene = pd.DataFrame(data= Dict( :genes=> plot_names, :c=>plot_cs, :cond=>plot_conds))
fig,ax = PyPlot.subplots(figsize=(15,7))
ax = sns.violinplot(x="genes", y="c", hue="cond",split=true, data=gene,palette="Set1" , inner = nothing)
ax[:legend](fontsize=16, loc= "upper right")
ax[:set_ylabel]("Product abunance",fontsize=16)
#here
ax[:set_xlabel]("New genes",fontsize=16)
ax[:tick_params](labelsize=16)
#here 2
name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/new_genes_real_data/new_genes_violine.pdf")
savefig(name)
