using PyPlot
using ABC_Bursting
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport pandas as pd
using PyCall, PyPlot; @pyimport numpy as np

# genes = ["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]
genes = [ "ENSMUSG00000012396"]
sampled_parameter_combinations = []
fileending = "_serum_0.1_fix_100_50000"
n_samples=100
#read experimental data
all_genes_min_max = []
all_genes_experimental_results = []
all_genes_all_simulation_results = []
all_genes_samples = []
data = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_serum_filtereddata.txt")
for name in genes
      pathtofile_lisi = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/500/serum/model1_",name,fileending,".txt")
      x=read_one_m1(pathtofile_lisi,  3, n_samples,3)
      print(x)
      cum_posterior = Vector(n_samples)
      cumsum!(cum_posterior, x[6])
      sample_rand = rand(500)
      winner_inds = []
      for s in sample_rand
          for i in 1:length(cum_posterior)
                if s<cum_posterior[i]
                  push!(winner_inds, [x[1][i], x[2][i], x[3][i]])
                  break
                end
          end
      end
      push!(sampled_parameter_combinations, winner_inds)
end
for i in 1:length(genes)
    gene = genes[i]
    #read and store experimental data for gene "gene"
    one_gene_experimental_result = Get_single_gene_RNAseq_data(gene,data)
    push!(all_genes_experimental_results,one_gene_experimental_result)
    #simulate and do histogram of samples
    one_gene_all_simulation_results= []
    one_gene_bins = Dict()
    one_gene_bins["bins"] =  range(1,20)
    for sim in 1:length(sampled_parameter_combinations[i])
        one_gene_one_simulation_result= generate_single_simulation_samples_m1_50000(sampled_parameter_combinations[i][sim], n_samples, data_scaling_RNAseq)
        push!(one_gene_all_simulation_results, one_gene_one_simulation_result)
    end
    min_ele = min(min(one_gene_all_simulation_results...)...)
    max_ele = max(max(one_gene_all_simulation_results...)...)
    push!(all_genes_min_max, [min_ele, max_ele])
    for  sim in 1:length(sampled_parameter_combinations[i])
        bins, indata =np.histogram(one_gene_all_simulation_results[sim],density=true,bins = 20, range=(min_ele,max_ele))
        one_gene_bins["sample_$(sim)"] = bins
    end
    push!(all_genes_samples,one_gene_bins)
    push!(all_genes_all_simulation_results,one_gene_all_simulation_results)
#     PyPlot.ioff()
#     fig, ax = PyPlot.subplots(figsize=(20,15),ncols=3, nrows=1)
#     sns.distplot(res_i,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "experimental",axlabel = "pruduct abundance",ax=ax[1,1])
#     sns.distplot(res_i,bins=20, kde= true, hist=false,rug= false ,  color = "orchid", label = "experimental",axlabel = "pruduct abundance",ax=ax[2,1])
#     for sim in 1:length(sampled_parameter_combinations[i])
#         sns.distplot(one_gene_one_simulation_result_samples[sim],bins=20, kde= true, hist=false,rug= false ,  color = "red", label = "simulated",ax=ax[2,1])
#     end
#     sns.distplot(one_gene_one_simulation_result_samples[1],bins=20, kde= true, hist=true,rug= true ,  color = "red", label = "simulated",axlabel = "pruduct abundance",ax=ax[3,1])
#     # name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/exper/model_1_cellcycle/", gene)
#     name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/testing-priors/obs-vs-sim/", gene,fileending,".png")
#     PyPlot.savefig(name)
end
for gene in 1:length(genes)
    data = all_genes_samples[gene]
    in_para =[]
    colnames = []
    for key in keys(all_genes_samples[gene])
        colname=(string(key))
        push!(colnames,colname)
    end
    df = pd.DataFrame(data, columns = colnames)
    for i in 1:length(colnames)
        push!(in_para,df[colnames[i]])
    end
    fig,ax = PyPlot.subplots()
      sns.tsplot(in_para,color="indianred", value="Frequency", time=np.linspace(all_genes_min_max[gene][1], all_genes_min_max[gene][2], 20),  ci=95,ax=ax)
      sns.distplot(all_genes_experimental_results[gene],norm_hist=true, bins=20, kde= true, hist=false,rug= false ,  color = "orchid", label = "experimental",axlabel = "Product abundance",ax=ax)

    name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/model-1-nanog/200/serum/", genes[gene],fileending,".png")
     savefig(name)
end
