using ABC_Bursting,PyCall, PyPlot
@pyimport seaborn as sns
@pyimport pandas as pd

gene_names = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/gene_names/genes_high_var_serum.txt",'\t','\n')
gene_names_en = gene_names[1:end,:,:]
real_names = Dict()
gene_names_en[2,3]
for gene_counter in 1:length(gene_names_en[:,1,1])
   real_names[gene_names_en[gene_counter]] = gene_names_en[gene_counter,2]
end

# genes = ["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]
# genes = [ "ENSMUSG00000041431"]
# colors = ["darkblue", "cornflowerblue", "orangered", "darkred"]

# real_names = ["Gli2 in serum","Nr0b1 in serum",  "Gli2 in 2i" , "Zfp42 in 2i" ]
# markers = ["+" , "*" , "^","+" , "*" , "^"]
sampled_parameter_combinations_2i = []
sampled_parameter_combinations_serum = []
n_samples = 500
n_tester = 100

fileending_2i="model_8_serum_genes_2i_data_"
fileending_serum="model_8_serum_genes_serum_data_"
gene ="ENSMUSG00000012396"
n_samples
path_to_file_2i = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_8_0.4_5/serum_genes/2i_data/",fileending_2i,gene,".txt")
path_to_file_serum = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_8_0.4_5/serum_genes/serum_data/",fileending_serum,gene,".txt")
data_file_2i=read_one_m8(path_to_file_2i,  3, n_samples,4)
data_file_serum=read_one_m8(path_to_file_serum,  3, n_samples,4)
#
# path_to_file_2i_m9 = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_9_huge/serum_genes/2i_data/model_9_serum_genes_2i_data_ENSMUSG00000012396_david.txt")
# path_to_file_serum_m9 = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_9_huge/serum_genes/serum_data/model_9_serum_genes_serum_data_ENSMUSG00000012396_david.txt")
# data_file_2i_m9=read_one_m9(path_to_file_2i_m9,  3, n_samples,4)
# data_file_serum_m9=read_one_m9(path_to_file_serum_m9,  3, n_samples,4)
#

PyPlot.ion()

      # sns.distplot(data_file_2i[4],bins=15, kde= true, hist=true,rug= false ,color = "salmon", label = "lisi")
      # sns.distplot(data_file_serum[4],bins=15, kde= true, hist=true,rug= false ,color = "cornflowerblue", label = "lisi")

            plot_cs = vcat(data_file_2i[1],data_file_serum[1],data_file_2i[2],data_file_serum[2],data_file_2i[3],data_file_serum[3],data_file_2i[4],data_file_serum[4])
            x1 = fill("2i",length(data_file_2i[4]))
            x2 = fill("serum",length(data_file_serum[4]))
            xes = vcat(x1,x2,x1,x2,x1,x2,x1,x2)
            aa=fill("Act. ",length(data_file_serum[4]))
            ab=fill("Deact.",length(data_file_serum[4]))
            ac=fill("Deg. ",length(data_file_serum[4]))
            ad=fill("Fb. k",length(data_file_serum[4]))
            genes = vcat(aa,aa,ab,ab,ac,ac,ad, ad)
            print(length(xes), length(genes),length(plot_cs))

            gene = pd.DataFrame(data= Dict( :Genes=> genes, :c=>(plot_cs), :Condition=>xes))
            fig,ax = PyPlot.subplots()
            ax = sns.boxplot(x="Genes", y="c",data=gene, hue="Condition", palette="Set1")
            ax[:tick_params](labelsize=16)
            ax[:legend](fontsize=16, loc ="upper left")
            ax[:yaxis][:set_label_position]("right")
            ax[:set_ylabel]("Parameter value",fontsize=16)
            ax[:set_xlabel]("",fontsize=16)
            name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/params/nanog_all_paramsr.pdf")

 PyPlot.savefig(name)

a
a=convert(Array{Float64,1},data_file_2i[4])
b=convert(Array{Float64,1},data_file_serum[4])
a2 =[10^i for i in a]
b2 =[10^i for i in b]
a = kolmogorov_smirnov_distance(a2,b2)

2^.([1,2])
