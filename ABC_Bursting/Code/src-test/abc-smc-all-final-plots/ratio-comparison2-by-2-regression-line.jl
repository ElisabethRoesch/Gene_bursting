using PyPlot
using ABC_Bursting
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport pandas as pd
# genes = ["ENSMUSG00000003032", "ENSMUSG00000012396", "ENSMUSG00000018604", "ENSMUSG00000020717", "ENSMUSG00000021255", "ENSMUSG00000021835", "ENSMUSG00000022528", "ENSMUSG00000025056", "ENSMUSG00000038793", "ENSMUSG00000048402", "ENSMUSG00000051176"]
# genes = ["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"] #2i genes
genes = [ "ENSMUSG00000003032"]
sampled_parameter_combinations_2i = []
sampled_parameter_combinations_serum = []
n_samples = 500
n_tester = 100

fileending_2i = "_0.1_500_250000_serum_genes_2i_data"
fileending_serum = "_0.1_500_250000_serum_genes_serum_data"


for gene_counter in 1:length(genes)
      name=genes[gene_counter]
      path_to_file_2i = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/explore_distr/model_1/all_high_var_in_serum_genes_in_2i/model1_",name,fileending_2i,".txt")
      path_to_file_serum = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/explore_distr/model_1/all_high_var_in_serum_genes_in_serum/model1_",name,fileending_serum,".txt")

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

      #plot


      acts_2i = []
      deacts_2i = []
      degs_2i = []
      for i in 1:length(sampled_parameter_combinations_2i[gene_counter])
            push!(acts_2i,sampled_parameter_combinations_2i[gene_counter][i][1])
            push!(deacts_2i,sampled_parameter_combinations_2i[gene_counter][i][2])
            push!(degs_2i,sampled_parameter_combinations_2i[gene_counter][i][3])
      end

      a=convert(Array{Float64,1},acts_2i)
      b=convert(Array{Float64,1},deacts_2i)
      c=convert(Array{Float64,1},degs_2i)
      a2= [10^n for n in a]
      b2= [10^n for n in b]
      c2= [10^n for n in c]
      r_deact_deg_2i =  log10.((b2)./(c2))
      r_act_deg_2i =  log10.((a2)./(c2))


      acts_serum = []
      deacts_serum = []
      degs_serum = []
      for i in 1:length(sampled_parameter_combinations_serum[gene_counter])
            push!(acts_serum,sampled_parameter_combinations_serum[gene_counter][i][1])
            push!(deacts_serum,sampled_parameter_combinations_serum[gene_counter][i][2])
            push!(degs_serum,sampled_parameter_combinations_serum[gene_counter][i][3])
      end
      a=convert(Array{Float64,1},acts_serum)
      b=convert(Array{Float64,1},deacts_serum)
      c=convert(Array{Float64,1},degs_serum)
      a2= [10^n for n in a]
      b2= [10^n for n in b]
      c2= [10^n for n in c]
      r_deact_deg_serum =  log10.((b2)./(c2))
      r_act_deg_serum =  log10.((a2)./(c2))
      centroid_x_serum ,centroid_y_serum = (sum(r_deact_deg_serum)/length(r_deact_deg_serum), sum(r_act_deg_serum)/length(r_act_deg_serum))
      centroid_x_2i ,centroid_y_2i = (sum(r_deact_deg_2i)/length(r_deact_deg_2i), sum(r_act_deg_2i)/length(r_act_deg_2i))
      println(      centroid_x_serum ,centroid_y_serum,centroid_x_2i ,centroid_y_2i)

      fig,ax = PyPlot.subplots()
            # PyPlot.scatter([-1,2.2],[0.,0.], color="white", marker = "*")
             PyPlot.scatter([centroid_x_serum],[centroid_y_serum], color= "cornflowerblue", marker = "o")
             PyPlot.scatter([centroid_x_2i],[centroid_y_2i], color="salmon", marker ="o")
            sns.regplot(r_deact_deg_serum,r_act_deg_serum ,color="#b1caf6",marker = "+",fit_reg =false)
            sns.regplot(r_deact_deg_2i,r_act_deg_2i, color = "salmon",marker = "+", fit_reg =false)
            ax[:set_ylabel]("log(act/deg)")
            ax[:set_xlabel]("log(deact/deg)")

      savefig(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/ratios-all-genes-base-model/",genes[gene_counter],fileending_serum,"ratios_center.png"))
end
rgba(255,0,0,0.3)
