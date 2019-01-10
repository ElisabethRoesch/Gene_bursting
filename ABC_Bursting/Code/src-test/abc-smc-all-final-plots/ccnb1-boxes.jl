using ABC_Bursting,PyCall, PyPlot
@pyimport seaborn as sns
@pyimport pandas as pd

n_samples=500
path_to_file_2i = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/cell_cycle/2i_genes/2i_data/model1_cell_cycle_ENSMUSG00000041431_2i_genes_2i_data_0.1_500_250000_16666_33332.txt")
path_to_file_serum = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/cell_cycle/2i_genes/serum_data/model1_cell_cycle_ENSMUSG00000041431_2i_genes_serum_data_0.1_500_250000_16666_33332.txt")
aa=1
data_file_2i=read_one_m1(path_to_file_2i,  3, n_samples,3)
data_file_serum=read_one_m1(path_to_file_serum,  3, n_samples,3)

PyPlot.ion()
            plot_cs = vcat(data_file_2i[1],data_file_serum[1],data_file_2i[2],data_file_serum[2],data_file_2i[3],data_file_serum[3])
            x1 = fill("2i",length(data_file_2i[1]))
            x2 = fill("serum",length(data_file_serum[1]))
            xes = vcat(x1,x2,x1,x2,x1,x2)
            aa=fill("Act.",length(data_file_serum[1]))
            ab=fill("Deact.",length(data_file_serum[1]))
            ac=fill("Deg.",length(data_file_serum[1]))
            genes = vcat(aa,aa,ab,ab,ac,ac)
            print(length(xes), length(genes),length(plot_cs))

            gene = pd.DataFrame(data= Dict( :Genes=> genes, :c=>(plot_cs), :Condition=>xes))
            fig,ax = PyPlot.subplots()
            ax = sns.boxplot(x="Genes", y="c",data=gene, hue="Condition", palette="Set1",whis=10.)
            # ax[:yaxis][:set_label_position]("right")
            ax[:tick_params](labelsize=16)

            ax[:legend](fontsize=16, loc ="upper right")
            ax[:set_ylabel]("Parameter value",fontsize=16)
            ax[:set_xlabel]("",fontsize=16)
            # name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/params/ccnb1_.pdf")

# PyPlot.savefig(name)
