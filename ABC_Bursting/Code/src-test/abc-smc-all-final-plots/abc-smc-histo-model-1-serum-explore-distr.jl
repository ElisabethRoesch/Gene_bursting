using PyPlot
using ABC_Bursting
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport pandas as pd
using PyCall, PyPlot; @pyimport numpy as np
# ENSMUSG00000025056 ENSMUSG00000048402 ENSMUSG00000021835 ENSMUSG00000048402
# genes = [ "ENSMUSG00000021835","ENSMUSG00000025056",
# genes =["ENSMUSG00000048402"]
genes = ["ENSMUSG00000006589", "ENSMUSG00000019539", "ENSMUSG00000020571", "ENSMUSG00000020954", "ENSMUSG00000022881", "ENSMUSG00000024646", "ENSMUSG00000025130", "ENSMUSG00000025359", "ENSMUSG00000047751", "ENSMUSG00000050953", "ENSMUSG00000062270"]

#in seurm ENSMUSG00000025056 ENSMUSG00000048402
# in 2i : ENSMUSG00000051176, ENSMUSG00000048402
# genes = ["ENSMUSG00000051176","ENSMUSG00000020717","ENSMUSG00000048402"]
sampled_parameter_combinations_serum = []
n_samples = 500
n_tester = 500
binsize = 50
#read experimental data_2i
all_genes_min_max_serum = []
all_genes_experimental_results_serum = []
all_genes_all_simulation_results_serum = []
all_genes_samples_serum = []
data_serum = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_serum_filtereddata.txt")
for name in genes
      path_to_file_serum =  string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/explore_distr/model_1/serum/model1_",name, "_serum_0.1_fix_500_250000.txt")
      data_file_serum=read_one_m1(path_to_file_serum,  3, n_samples,3)
      cum_posterior_serum = Vector(n_samples)
      cumsum!(cum_posterior_serum, data_file_serum[6])
      sample_rand_serum = rand(n_tester)
      winner_inds_serum = []
      for s in sample_rand_serum
          for i in 1:length(cum_posterior_serum)
                if s<cum_posterior_serum[i]
                  push!(winner_inds_serum, [data_file_serum[1][i], data_file_serum[2][i], data_file_serum[3][i]])
                  break
                end
          end
      end
      push!(sampled_parameter_combinations_serum, winner_inds_serum)
end
for i in 1:length(genes)
    gene = genes[i]
    #read and store experimental data_2i for gene "gene"
    one_gene_experimental_result_serum = Get_single_gene_RNAseq_data(gene,data_serum)
    push!(all_genes_experimental_results_serum,one_gene_experimental_result_serum)
    #simulate and do histogram of samples
    one_gene_all_simulation_results_serum= []
    one_gene_bins_serum = Dict()
    one_gene_bins_serum["bins"] =  range(1,binsize)
    for sim in 1:length(sampled_parameter_combinations_serum[i])
        one_gene_one_simulation_result_serum= generate_single_simulation_samples_m1_50000(sampled_parameter_combinations_serum[i][sim], n_samples, data_scaling_RNAseq)
        push!(one_gene_all_simulation_results_serum, one_gene_one_simulation_result_serum)
    end
    min_ele_serum = min(min(one_gene_all_simulation_results_serum...)...)
    max_ele_serum = max(max(one_gene_all_simulation_results_serum...)...)
    push!(all_genes_min_max_serum, [min_ele_serum, max_ele_serum])
    for  sim in 1:length(sampled_parameter_combinations_serum[i])
        bins_serum, indata_serum =np.histogram(one_gene_all_simulation_results_serum[sim], density = true, bins = binsize, range = (0, max_ele_serum))
        one_gene_bins_serum["sample_$(sim)"] = bins_serum
    end
    push!(all_genes_samples_serum,one_gene_bins_serum)
    push!(all_genes_all_simulation_results_serum,one_gene_all_simulation_results_serum)
end

for gene in 1:length(genes)
    data_serum = all_genes_samples_serum[gene]

    in_para_serum =[]

    colnames_serum = []

    for key in keys(all_genes_samples_serum[gene])
        colname=(string(key))
        push!(colnames_serum,colname)
    end

    df_serum = pd.DataFrame(data_serum, columns = colnames_serum)


    for i in 1:length(colnames_serum)
        push!(in_para_serum,df_serum[colnames_serum[i]])
    end

    fig,ax = PyPlot.subplots()

      sns.tsplot(in_para_serum,color="darkblue", time=np.linspace(0, all_genes_min_max_serum[gene][2], binsize),  ci=95,ax=ax)
      sns.distplot(all_genes_experimental_results_serum[gene],norm_hist=true, bins=binsize, kde= false, hist=true,rug= false ,  color = "grey", label = "Experiment")
      sns.distplot([0.],norm_hist=true, bins=binsize, kde= false, hist=false,rug= true ,  color = "darkblue", label = "Simulation")

      ax[:set_ylabel]("Frequency",fontsize=16)
      ax[:set_xlabel]("Product abundance",fontsize=16)
      ax[:legend](fontsize=16,loc="upper right")
      name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/explore-distr/final/", genes[gene],"_n_tester",string(n_tester),"_binsize",string(binsize),"serum.pdf")
     savefig(name)
end
