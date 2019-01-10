using ABC_Bursting
using PyCall, PyPlot; @pyimport seaborn as sns
# genes = ["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]
# genes = ["ENSMUSG00000027715", "ENSMUSG00000029177"]
genes = [ "ENSMUSG00000006398"]
genes=["ENSMUSG00000006398","ENSMUSG00000027715"]
#, "ENSMUSG00000029177"]

#sample parameter combinations of ABC results
xes=[]
winners=[]
fileending= "_2i_0.1_fix_100_50000"

for name in genes
    #pathtofile_lisi = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_1_cellcycle_2i_medium/model_1_cellcycle_medium_",name,"_2i_0.8_100.txt")

  pathtofile_lisi = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/model1_",name,fileending,".txt")
  x=read_one_m1(pathtofile_lisi,  3, 100,3)
  push!(xes,x)
  print(x)
  cum_posterior = Vector(100)
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
distrs= []
data = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_2i_filtereddata.txt")
for i in 1:length(genes)
    gene = genes[i]
    res_i = Get_single_gene_RNAseq_data(gene,data)
    push!(res,res_i)
    resb = []
    for sim in 1:length(winners[i])
        # resb_sim = cell_cycle_simulation_generate_single_simulation_samples_m1(winners[i][sim], 100, data_scaling_RNAseq)
        resb_sim = generate_single_simulation_samples_m1_200000(winners[i][sim], 200, data_scaling_RNAseq)
        push!(resb,resb_sim)
    end
    push!(res_sim,resb)
    PyPlot.ioff()
    fig, ax = PyPlot.subplots(figsize=(20,15),ncols=3, nrows=1)
    sns.distplot(res_i,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "experimental",axlabel = "pruduct abundance",ax=ax[1,1])
    sns.distplot(res_i,bins=20, kde= true, hist=false,rug= false ,  color = "orchid", label = "experimental",axlabel = "pruduct abundance",ax=ax[2,1])
    for sim in 1:length(winners[i])
        sns.distplot(resb[sim],bins=20, kde= true, hist=false,rug= false ,  color = "red", label = "simulated",ax=ax[2,1])
    end
    sns.distplot(resb[1],bins=20, kde= true, hist=true,rug= true ,  color = "red", label = "simulated",axlabel = "pruduct abundance",ax=ax[3,1])
    # name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/exper/model_1_cellcycle/", gene)
    name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/testing-priors/obs-vs-sim/", gene,fileending,".png")

    PyPlot.savefig(name)
end

#
# for name in genes
#   pathtofile_lisi = string(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_8_2i_medium/model_8_easy_",name),"_2i_0.8_100.txt")
#   x=read_one_m8(pathtofile_lisi, pathtofile_lisi, 3, 100, 5)
#   push!(distrs,x)
# end
#
