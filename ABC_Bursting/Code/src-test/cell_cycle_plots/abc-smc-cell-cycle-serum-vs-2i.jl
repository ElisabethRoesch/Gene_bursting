using PyPlot
using ABC_Bursting
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport pandas as pd
using PyCall, PyPlot; @pyimport numpy as np

# genes = ["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]
genes = [ "ENSMUSG00000027715", "ENSMUSG00000041431"]

sampled_parameter_combinations_2i = []
sampled_parameter_combinations_serum = []

n_samples = 500
n_tester = 500
tau = 500
tmax = tau*n_samples
#TODO here twice
t1 = 16666
t2 =  33332
fileending_2i = string("_2i_genes_2i_data_0.1_500_250000_",t1,"_",t2)
fileending_serum = string("_2i_genes_serum_data_0.1_500_250000_",t1,"_",t2)
binsize = 15
#read experimental data_2i
all_genes_min_max_2i = []
all_genes_min_max_serum = []
all_genes_min_max_all = []


all_genes_experimental_results_2i = []
all_genes_experimental_results_serum = []

all_genes_all_simulation_results_2i = []
all_genes_all_simulation_results_serum = []

all_genes_samples_2i = []
all_genes_samples_serum = []


data_2i = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_2i_filtereddata.txt")
data_serum = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_serum_filtereddata.txt")

for name in genes
      path_to_file_2i = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/cell_cycle/2i_genes/2i_data/model1_cell_cycle_",name,fileending_2i,".txt")
      path_to_file_serum = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/cell_cycle/2i_genes/serum_data/model1_cell_cycle_",name,fileending_serum,".txt")

      data_file_2i=read_one_m1(path_to_file_2i,  3, n_samples,3)
      data_file_serum=read_one_m1(path_to_file_serum,  3, n_samples,3)

      cum_posterior_2i = Vector(n_samples)
      cum_posterior_serum = Vector(n_samples)

      cumsum!(cum_posterior_2i, data_file_2i[6])
      cumsum!(cum_posterior_serum, data_file_serum[6])

      sample_rand_2i = rand(n_tester)
      sample_rand_serum = rand(n_tester)

      winner_inds_2i = []
      winner_inds_serum = []

      for s in sample_rand_2i
          for i in 1:length(cum_posterior_2i)
                if s<cum_posterior_2i[i]
                  push!(winner_inds_2i, [data_file_2i[1][i], data_file_2i[2][i], data_file_2i[3][i]])
                  break
                end
          end
      end

      for s in sample_rand_serum
          for i in 1:length(cum_posterior_serum)
                if s<cum_posterior_serum[i]
                  push!(winner_inds_serum, [data_file_serum[1][i], data_file_serum[2][i], data_file_serum[3][i]])
                  break
                end
          end
      end

      push!(sampled_parameter_combinations_2i, winner_inds_2i)
      push!(sampled_parameter_combinations_serum, winner_inds_serum)

end
for i in 1:length(genes)
    gene = genes[i]
    #read and store experimental data_2i for gene "gene"
    one_gene_experimental_result_2i = Get_single_gene_RNAseq_data(gene,data_2i)
    one_gene_experimental_result_serum = Get_single_gene_RNAseq_data(gene,data_serum)

    push!(all_genes_experimental_results_2i,one_gene_experimental_result_2i)
    push!(all_genes_experimental_results_serum,one_gene_experimental_result_serum)

    #simulate and do histogram of samples
    one_gene_all_simulation_results_2i= []
    one_gene_all_simulation_results_serum= []

    one_gene_bins_2i = Dict()
    one_gene_bins_serum = Dict()

    one_gene_bins_2i["bins"] =  range(1,binsize)
    one_gene_bins_serum["bins"] =  range(1,binsize)

    for sim in 1:length(sampled_parameter_combinations_2i[i])
        one_gene_one_simulation_result_2i = cell_cycle_simulation_generate_single_simulation_samples_m1_16666_33332(sampled_parameter_combinations_2i[i][sim], n_samples, data_scaling_RNAseq)
        push!(one_gene_all_simulation_results_2i, one_gene_one_simulation_result_2i)
    end
    for sim in 1:length(sampled_parameter_combinations_serum[i])
        one_gene_one_simulation_result_serum= cell_cycle_simulation_generate_single_simulation_samples_m1_16666_33332(sampled_parameter_combinations_serum[i][sim], n_samples, data_scaling_RNAseq)
        push!(one_gene_all_simulation_results_serum, one_gene_one_simulation_result_serum)
    end

    min_ele_2i = min(min(one_gene_all_simulation_results_2i...)...)
    max_ele_2i = max(max(one_gene_all_simulation_results_2i...)...)

    min_ele_serum = min(min(one_gene_all_simulation_results_serum...)...)
    max_ele_serum = max(max(one_gene_all_simulation_results_serum...)...)
    min_ele_all = min(min_ele_2i,min_ele_serum)
    max_ele_all = max(max_ele_2i,max_ele_serum)
    push!(all_genes_min_max_2i, [min_ele_2i, max_ele_2i])
    push!(all_genes_min_max_serum, [min_ele_serum, max_ele_serum])

    push!(all_genes_min_max_all,[min_ele_all, max_ele_all])

    for  sim in 1:length(sampled_parameter_combinations_2i[i])
        bins_2i, indata_2i =np.histogram(one_gene_all_simulation_results_2i[sim],density=true,bins = binsize, range=(min_ele_all,max_ele_all))
        one_gene_bins_2i["sample_$(sim)"] = bins_2i
    end

    for  sim in 1:length(sampled_parameter_combinations_serum[i])
        bins_serum, indata_serum =np.histogram(one_gene_all_simulation_results_serum[sim],density=true,bins = binsize, range=(min_ele_all,max_ele_all))
        one_gene_bins_serum["sample_$(sim)"] = bins_serum
    end

    push!(all_genes_samples_2i,one_gene_bins_2i)
    push!(all_genes_all_simulation_results_2i,one_gene_all_simulation_results_2i)

    push!(all_genes_samples_serum,one_gene_bins_serum)
    push!(all_genes_all_simulation_results_serum,one_gene_all_simulation_results_serum)
end

for gene in 1:length(genes)
    data_2i = all_genes_samples_2i[gene]
    data_serum = all_genes_samples_serum[gene]

    in_para_2i =[]
    in_para_serum =[]

    colnames_2i = []
    colnames_serum = []

    for key in keys(all_genes_samples_2i[gene])
        colname=(string(key))
        push!(colnames_2i,colname)
    end

    for key in keys(all_genes_samples_serum[gene])
        colname=(string(key))
        push!(colnames_serum,colname)
    end

    df_2i = pd.DataFrame(data_2i, columns = colnames_2i)
    df_serum = pd.DataFrame(data_serum, columns = colnames_serum)

    for i in 1:length(colnames_2i)
        push!(in_para_2i,df_2i[colnames_2i[i]])
    end

    for i in 1:length(colnames_serum)
        push!(in_para_serum,df_serum[colnames_serum[i]])
    end

    fig,ax = PyPlot.subplots()

      sns.tsplot(in_para_serum,color="colorblind", value="Frequency", time=np.linspace(all_genes_min_max_all[gene][1], all_genes_min_max_all[gene][2], binsize),  ci=95,ax=ax)
      sns.tsplot(in_para_2i,color="salmon", value="Frequency", time=np.linspace(all_genes_min_max_all[gene][1], all_genes_min_max_all[gene][2], binsize),  ci=95,ax=ax)
      sns.distplot(all_genes_experimental_results_2i[gene],norm_hist=true, bins=binsize, kde= true, hist=false,rug= false ,  color = "red", label = "experimental 2i",axlabel = "Product abundance",ax=ax)
      sns.distplot(all_genes_experimental_results_serum[gene],norm_hist=true, bins=binsize, kde= true, hist=false,rug= false ,  color = "blue", label = "experimental serum",axlabel = "Product abundance",ax=ax)
    name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/posterior_quantile/cc/", genes[gene],fileending_2i,fileending_serum,"_n_tester",string(n_tester),"_binsize",string(binsize),"_t1_",string(t1),"_t2_",string(t2),".png")
     savefig(name)
end
