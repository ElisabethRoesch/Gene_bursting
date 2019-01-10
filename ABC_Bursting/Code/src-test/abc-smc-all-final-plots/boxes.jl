using ABC_Bursting, PyCall, PyPlot
@pyimport seaborn as sns
@pyimport pandas as pd
# @pyimport scipy as spyy
# "ENSMUSG00000071547" not found
@pyimport scipy.stats as sta
# genes = ["ENSMUSG00000003032", "ENSMUSG00000012396", "ENSMUSG00000018604", "ENSMUSG00000020717", "ENSMUSG00000021255", "ENSMUSG00000021835", "ENSMUSG00000022528", "ENSMUSG00000025056", "ENSMUSG00000038793", "ENSMUSG00000048402", "ENSMUSG00000051176"]
# genes = ["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]
genes = ["ENSMUSG00000006589", "ENSMUSG00000019539", "ENSMUSG00000020571", "ENSMUSG00000020954", "ENSMUSG00000022881", "ENSMUSG00000024646", "ENSMUSG00000025130", "ENSMUSG00000025359", "ENSMUSG00000047751", "ENSMUSG00000050953", "ENSMUSG00000062270"]
all_genes_length = length(genes)
gene_names = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/gene_names/new_genes.txt",'\t','\n')
gene_names_en = gene_names[1:end,:,:]
real_names = Dict()
gene_names_en[2,3]
for gene_counter in 1:length(gene_names_en[:,1,1])
   real_names[gene_names_en[gene_counter]] = gene_names_en[gene_counter,2]
end

# genes = [ "ENSMUSG00000041431"]
# colors = ["darkblue", "cornflowerblue", "orangered", "darkred"]

# real_names = ["Gli2 in serum","Nr0b1 in serum",  "Gli2 in 2i" , "Zfp42 in 2i" ]
# markers = ["+" , "*" , "^","+" , "*" , "^"]
sampled_parameter_combinations_2i = []
sampled_parameter_combinations_serum = []
n_samples = 100
n_tester = 500

acts = []
deacts = []
degs = []
acts_2i = []
deacts_2i = []
degs_2i = []
# endfile = ["serum","serum","2i","2i"]
counter = 0
PyPlot.ioff()

for name in genes
      counter=counter+1
      path_to_file_2i =  string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model-1-new-genes/2i/model1_2i_",name, "_new_gene.txt")
      path_to_file_serum =  string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model-1-new-genes/serum/model1_serum_", name, "_new_gene.txt.txt")
      data_file_2i=read_one_m1(path_to_file_2i,  3, n_samples,3)
      data_file_serum=read_one_m1(path_to_file_serum,  3, n_samples,3)
      cum_posterior_2i = Vector(n_samples)
      cum_posterior_serum = Vector(n_samples)
      cumsum!(cum_posterior_serum, data_file_serum[6])
      cumsum!(cum_posterior_2i, data_file_2i[6])
      sample_rand_serum = rand(n_tester)
      sample_rand_2i = rand(n_tester)
      winner_inds_serum = []
      winner_inds_2i = []
      for s in sample_rand_serum
          for i in 1:length(cum_posterior_serum)
                if s<cum_posterior_serum[i]
                  push!(winner_inds_serum, [data_file_serum[1][i], data_file_serum[2][i], data_file_serum[3][i]])
                  break
                end
          end
      end
      for s in sample_rand_2i
          for i in 1:length(cum_posterior_2i)
                if s<cum_posterior_2i[i]
                  push!(winner_inds_2i, [data_file_2i[1][i], data_file_2i[2][i], data_file_2i[3][i]])
                  break
                end
          end
      end
      push!(sampled_parameter_combinations_2i, winner_inds_2i)

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
      push!(acts, a2)
      push!(deacts, b2)
      push!(degs, c2)
end

for name in 1:length(genes)
      actst = []
      deactst = []
      degst = []
      for i in 1:length(sampled_parameter_combinations_2i[name])
            push!(actst,convert(Float64,sampled_parameter_combinations_2i[name][i][1]))
            push!(deactst,convert(Float64,sampled_parameter_combinations_2i[name][i][2]))
            push!(degst,convert(Float64,sampled_parameter_combinations_2i[name][i][3]))
      end
      a=convert(Array{Float64,1},actst)
      b=convert(Array{Float64,1},deactst)
      c=convert(Array{Float64,1},degst)
      a2= [10^n for n in a]
      b2= [10^n for n in b]
      c2= [10^n for n in c]
      push!(acts_2i, a2)
      push!(deacts_2i, b2)
      push!(degs_2i, c2)
end

aks=[]
deaks=[]
degaks=[]
pyaks=[]
pydeaks=[]
pydegaks=[]

countergene = 0
for name in genes
      countergene=countergene+1
      push!(pyaks,sta.ks_2samp(acts[countergene],acts_2i[countergene]))
      push!(pydeaks,sta.ks_2samp(deacts[countergene],deacts_2i[countergene]))
      push!(pydegaks,sta.ks_2samp(degs[countergene],degs_2i[countergene]))

      push!(aks,kolmogorov_smirnov_distance(acts[countergene],acts_2i[countergene]))
      push!(deaks,kolmogorov_smirnov_distance(deacts[countergene],deacts_2i[countergene]))
      push!(degaks,kolmogorov_smirnov_distance(degs[countergene],degs_2i[countergene]))
      plot_cs = vcat(acts[countergene],deacts[countergene],degs[countergene],acts_2i[countergene],deacts_2i[countergene],degs_2i[countergene])
      x1 = fill("Activation rate",length(acts[countergene]))
      x2 = fill("Deactivation rate",length(acts[countergene]))
      x3 = fill("Degradation rate",length(acts[countergene]))
      xes = vcat(x1,x2,x3,x1,x2,x3)
      c1 = fill("serum",3*length(acts[countergene]))
      c2 = fill("2i",3*length(acts[countergene]))
      condi = vcat(c1,c2)


      gene = pd.DataFrame(data= Dict( :Parameters=> xes, :c=>log10.(plot_cs), :Condition=>condi))
      fig,ax = PyPlot.subplots(figsize=(15,10))
      ax = sns.boxplot(x="Condition", y="c",data=gene, hue="Parameters")
      ax[:legend](fontsize=16)
      ax[:set_ylabel]("Product abunance",fontsize=16)
      ax[:set_xlabel]("",fontsize=16)
      my_label=string(real_names[genes[countergene]])

      #here
      ax[:tick_params](labelsize=16)
      #here 2
      name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/params/new_genes/",my_label,"_boxplot_sorted_cond.pdf")
      savefig(name)

end
println(pyaks,"\n",pydeaks,"\n",pydegaks)
println(aks,"\n",deaks,"\n",degaks)


#only one serum
# countergene = 0
# for name in genes
#       countergene=countergene+1
#       plot_cs = vcat(acts[countergene],deacts[countergene],degs[countergene])
#       x1 = fill("Activation rate",length(acts[countergene]))
#       x2 = fill("Deactivation rate",length(acts[countergene]))
#       x3 = fill("Degradation rate",length(acts[countergene]))
#       xes = vcat(x1,x2,x3)
#       c1 = fill("Serum condition of",3*length(acts[countergene]))
#       # c2 = fill("2i",3*length(acts[countergene]))
#       condi = c1
#       gene = pd.DataFrame(data= Dict( :Parameters=> xes, :c=>log10.(plot_cs), :Condition=>condi))
#       fig,ax = PyPlot.subplots(figsize=(10,10))
#       ax = sns.boxplot(x="Condition", y="c",data=gene, hue="Parameters")
#       ax[:legend](fontsize=16)
#       ax[:set_ylabel]("Product abunance",fontsize=16)
#       my_label=string(real_names[genes[countergene]])
#       ax[:set_xlabel](my_label,fontsize=16)
#       #here
#       ax[:tick_params](labelsize=16)
#       #here 2
#       name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/params/",my_label,"_logged.pdf")
#       savefig(name)
# end
# #
# countergene = 0
# for name in genes
#       countergene=countergene+1
#       plot_cs = vcat(acts_2i[countergene],deacts_2i[countergene],degs_2i[countergene])
#       x1 = fill("Activation rate",length(acts_2i[countergene]))
#       x2 = fill("Deactivation rate",length(acts_2i[countergene]))
#       x3 = fill("Degradation rate",length(acts_2i[countergene]))
#       xes = vcat(x1,x2,x3)
#       c1 = fill("2i condition of",3*length(acts_2i[countergene]))
#       # c2 = fill("2i",3*length(acts[countergene]))
#       condi = c1
#
#
#       gene = pd.DataFrame(data= Dict( :Parameters=> xes, :c=>log10.(plot_cs), :Condition=>condi))
#       fig,ax = PyPlot.subplots(figsize=(10,10))
#       ax = sns.boxplot(x="Condition", y="c",data=gene, hue="Parameters")
#       ax[:legend](fontsize=16)
#       ax[:set_ylabel]("Product abunance",fontsize=16)
#       my_label=string(real_names[genes[countergene]])
#       ax[:set_xlabel](my_label,fontsize=16)
#
#       #here
#       ax[:tick_params](labelsize=16)
#       #here 2
#       print(kolmogorov_smirnov_distance(a2,b2))
#       name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/params/",my_label,"_logged_2i.pdf")
#       # savefig(name)
#
# end
