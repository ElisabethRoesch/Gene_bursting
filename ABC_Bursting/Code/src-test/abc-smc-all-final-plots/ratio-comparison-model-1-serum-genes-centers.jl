using PyPlot
using ABC_Bursting
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport pandas as pd
@pyimport matplotlib.patches as patches

#2i "ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431",
#serum genes = ["ENSMUSG00000003032", "ENSMUSG00000012396", "ENSMUSG00000018604", "ENSMUSG00000020717", "ENSMUSG00000021255", "ENSMUSG00000021835", "ENSMUSG00000022528", "ENSMUSG00000025056", "ENSMUSG00000038793", "ENSMUSG00000048402", "ENSMUSG00000051176"]

# genes = ["ENSMUSG00000003032", "ENSMUSG00000012396", "ENSMUSG00000018604", "ENSMUSG00000020717", "ENSMUSG00000021255", "ENSMUSG00000021835", "ENSMUSG00000022528", "ENSMUSG00000025056", "ENSMUSG00000038793", "ENSMUSG00000048402", "ENSMUSG00000051176"]
xes = ["ENSMUSG00000003032", "ENSMUSG00000012396", "ENSMUSG00000018604", "ENSMUSG00000020717", "ENSMUSG00000021255", "ENSMUSG00000022528", "ENSMUSG00000025056", "ENSMUSG00000038793", "ENSMUSG00000048402", "ENSMUSG00000051176"]
genes =vcat(xes,xes)
all_genes_length = length(xes)

gene_names = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/gene_names/genes_high_var_serum.txt",'\t','\n')
gene_names_en = gene_names[1:end,:,:]
real_names = Dict()
gene_names_en[2,3]
for gene_counter in 1:length(gene_names_en[:,1,1])
   real_names[gene_names_en[gene_counter]] = gene_names_en[gene_counter,2]
end

sampled_parameter_combinations_serum = []
n_samples = 500
n_tester = 500

all_r_deact_deg_serum = []
all_r_act_deg_serumm = []
all_centroid_x_serum = []
all_centroid_y_serum = []

genekind = ["serum","serum","serum","serum","serum","serum","serum","serum","serum","serum"]
genekind =vcat(genekind,genekind)
datakind2 = ["serum","serum","serum","serum","serum","serum","serum","serum","serum","serum"]
datakind = ["2i","2i","2i","2i","2i","2i","2i","2i","2i","2i"]
datakind=vcat(datakind,datakind2)
colors=[]
markers=[]
all_my_labels=[]
for i in 1:length(genekind)
      if genekind[i] == "2i"
            push!(markers,"o")
            if datakind[i] =="2i"
                  push!(colors,"salmon")
            else
                  push!(colors,"cornflowerblue")
            end
      else
            push!(markers,"+")
            if datakind[i] =="2i"
                  push!(colors,"salmon")
            else
                  push!(colors,"cornflowerblue")
            end
      end
end
counter = 0
for name in genes
      counter=counter+1
      path_to_file_serum =  string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/explore_distr/model_1/all_high_var_in_",genekind[counter],"_genes_in_",datakind[counter],"/model1_",name,"_0.1_500_250000_",genekind[counter],"_genes_",datakind[counter],"_data.txt")
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
for name in 1:length(genes)
      acts_serum = []
      deacts_serum = []
      degs_serum = []
      for i in 1:length(sampled_parameter_combinations_serum[name])
            push!(acts_serum,convert(Float64,sampled_parameter_combinations_serum[name][i][1]))
            push!(deacts_serum,convert(Float64,sampled_parameter_combinations_serum[name][i][2]))
            push!(degs_serum,convert(Float64,sampled_parameter_combinations_serum[name][i][3]))
      end
      a=convert(Array{Float64,1},acts_serum)
      b=convert(Array{Float64,1},deacts_serum)
      c=convert(Array{Float64,1},degs_serum)
      a2= [10^n for n in a]
      b2= [10^n for n in b]
      c2= [10^n for n in c]
      r_deact_deg_serum =  log10.((b2)./(c2))
      r_act_deg_serum =  log10.((a2)./(c2))
      push!(all_r_deact_deg_serum, r_deact_deg_serum)
      push!(all_r_act_deg_serumm, r_act_deg_serum)
      centroid_x_serum ,centroid_y_serum = (sum(r_deact_deg_serum)/length(r_deact_deg_serum), sum(r_act_deg_serum)/length(r_act_deg_serum))
      push!(all_centroid_x_serum, centroid_x_serum)
      push!(all_centroid_y_serum, centroid_y_serum)
      my_label=string(real_names[genes[name]])
      push!(all_my_labels, my_label)
end
fig,ax = PyPlot.subplots(figsize=(15,10))
      for name_counter in 1:all_genes_length
            PyPlot.plot([all_centroid_x_serum[name_counter],all_centroid_x_serum[name_counter+all_genes_length]],[all_centroid_y_serum[name_counter],all_centroid_y_serum[name_counter+all_genes_length]],"-",color="grey")
            ax[:annotate](all_my_labels[name_counter],xy=(all_centroid_x_serum[name_counter], all_centroid_y_serum[name_counter]), xytext=(all_centroid_x_serum[name_counter]+0.02,all_centroid_y_serum[name_counter]+0.02), fontsize = 15,zorder=100)
      end
      leftest = min(all_centroid_x_serum[1:all_genes_length]...)
      rightest = max(all_centroid_x_serum[1:all_genes_length]...)
      lowest  =  min(all_centroid_y_serum[1:all_genes_length]...)
      higherst = max(all_centroid_y_serum[1:all_genes_length]...)
      length_hori = rightest - leftest
      length_verti = higherst - lowest
      center_x = leftest + length_hori/(2)
      center_y = lowest + length_verti/(2)

      leftest_serum = min(all_centroid_x_serum[all_genes_length+1:2*all_genes_length]...)
      rightest_serum = max(all_centroid_x_serum[all_genes_length+1:2*all_genes_length]...)
      lowest_serum  =  min(all_centroid_y_serum[all_genes_length+1:2*all_genes_length]...)
      higherst_serum = max(all_centroid_y_serum[all_genes_length+1:2*all_genes_length]...)
      length_hori_serum = rightest_serum - leftest_serum
      length_verti_serum = higherst_serum - lowest_serum
      center_x_serum = leftest_serum + length_hori_serum/(2)
      center_y_serum = lowest_serum + length_verti_serum/(2)

      print("serum",(center_x_serum, center_y_serum), length_hori_serum+0.2, length_verti_serum+0.2 )
      print("\n2i",(center_x, center_y), length_hori+0.2, length_verti+0.2)
      # println("length_hori: ",length_hori, "length_verti: ",length_verti, "x: ",center_x, "y: ",center_y)
      cloud_serum = patches.Ellipse((center_x_serum, center_y_serum), length_hori_serum+0.2, length_verti_serum+0.2, facecolor="b", alpha=0.1)
      cloud_2i = patches.Ellipse((center_x, center_y), length_hori+0.2, length_verti+0.2, facecolor="r", alpha=0.1)
      ax[:add_artist](cloud_serum)
      ax[:add_artist](cloud_2i)
      PyPlot.plot(all_centroid_x_serum[1:all_genes_length],all_centroid_y_serum[1:all_genes_length], label="2i",color = "salmon",marker = "o",markersize=15, linestyle="")
      PyPlot.plot(all_centroid_x_serum[all_genes_length+1:all_genes_length*2],all_centroid_y_serum[all_genes_length+1:all_genes_length*2], label="serum",color = "cornflowerblue",marker = "o",markersize=15, linestyle="")
      ax[:legend](fontsize=16, loc ="lower right")
      ax[:set_ylabel]("log(act/deg)",fontsize=16)
      ax[:set_xlabel]("log(deact/deg)",fontsize=16)
      ax[:set_xlabel]("log(deact/deg)",fontsize=16)
      ax[:tick_params](labelsize=16)
      ax[:set_ylim]([-1.5,0.2])

      ax[:set_xlim]([-0.2,1.5])

# savefig(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/validate/serum_genes_ratios.pdf"))
