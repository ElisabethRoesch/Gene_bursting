using ABC_Bursting
using PyCall, PyPlot; @pyimport seaborn as sns

name = "ENSMUSG00000012396"
locations = [
            # "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/500/2i/model1_ENSMUSG00000012396_2i_0.1_fix_100_50000.txt",
            # "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/500/2i/model1_ENSMUSG00000012396_2i_0.1_fix_200_100000.txt",
            "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/500/2i/model1_ENSMUSG00000012396_2i_0.1_fix_500_250000.txt",
            # "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/500/serum/model1_ENSMUSG00000012396_2i_0.1_fix_100_50000.txt",
            # "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/500/serum/model1_ENSMUSG00000012396_2i_0.1_fix_200_100000.txt",
            # "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/500/serum/model1_ENSMUSG00000012396_2i_0.1_fix_500_250000.txt",
            # "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/1000/2i/model1_ENSMUSG00000012396_2i_0.1_fix_100_100000.txt",
            # "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/1000/2i/model1_ENSMUSG00000012396_2i_0.1_fix_200_200000.txt",
            "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/1000/2i/model1_ENSMUSG00000012396_2i_0.1_fix_500_500000.txt"]
            # "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/1000/serum/model1_ENSMUSG00000012396_2i_0.1_fix_100_100000.txt",
            # "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/1000/serum/model1_ENSMUSG00000012396_2i_0.1_fix_200_200000.txt",
            # "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/1000/serum/model1_ENSMUSG00000012396_2i_0.1_fix_500_500000.txt"]
# n_samples = []
locate_plots = [
            "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/model-1-nanog/500/2i/model1_ENSMUSG00000012396_2i_0.1_fix_500_250000.png",
            "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/model-1-nanog/1000/2i/model1_ENSMUSG00000012396_2i_0.1_fix_500_500000.png"
                ]
#sample parameter combinations of ABC results
xes=[]
winners=[]
for location in locations
    #pathtofile_lisi = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_1_cellcycle_2i_medium/model_1_cellcycle_medium_",name,"_2i_0.8_100.txt")

  pathtofile_lisi = location
  x=read_one_m1(pathtofile_lisi,  3, 500, 4)
  push!(xes,x)
  cum_posterior = Vector(500)
  cumsum!(cum_posterior, x[6])
  sample_rand = rand(10)
  winner_inds = []
  for s in sample_rand
      for i in 1:length(cum_posterior)
            if s<cum_posterior[i]
              push!(winner_inds, [x[1][i], x[2][i], x[3][i]])
              break
            end
      end
  end
  push!(winners,winner_inds)
end


#read experimental data
res = []
res_sim = []

data = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_2i_filtereddata.txt")
for i in 1:length(locations)

    res_i = Get_single_gene_RNAseq_data(name,data)
    push!(res,res_i)
    resb = []
    for sim in 1:length(winners[i])
        # resb_sim = cell_cycle_simulation_generate_single_simulation_samples_m1(winners[i][sim], 100, data_scaling_RNAseq)
        resb_sim = generate_single_simulation_samples_m1(winners[i][sim], 500, data_scaling_RNAseq)
        push!(resb,resb_sim)
    end
    push!(res_sim,resb)

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
      sns.tsplot(in_para,color="indianred", value="frequency", time=np.linspace(all_genes_min_max[gene][1], all_genes_min_max[gene][2], 20),  ci=95,ax=ax)
      sns.distplot(all_genes_experimental_results[gene],norm_hist=true, bins=20, kde= true, hist=false,rug= false ,  color = "orchid", label = "experimental",axlabel = "pruduct abundance",ax=ax)

    name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/posterior_quantile/", genes[gene],fileending,".png")
     savefig("ha")
end
