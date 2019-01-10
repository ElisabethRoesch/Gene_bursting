using PyPlot
using ABC_Bursting
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport pandas as pd

# genes = ["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]
genes = [ "ENSMUSG00000012396"]
sampled_parameter_combinations_2i = []
sampled_parameter_combinations_serum = []
n_samples = 500
n_tester = 500
tau = 500
tmax= tau*n_samples
fileending_2i = string("_2i_0.1_fix_",n_samples,"_",tmax)
fileending_serum = string("_serum_0.1_fix_",n_samples,"_",tmax)





for name in genes
      path_to_file_2i = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/",tau,"/2i/model1_",name,fileending_2i,".txt")
      path_to_file_serum = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/",tau,"/serum/model1_",name,fileending_serum,".txt")

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



acts_2i = []
deacts_2i = []
degs_2i = []
for i in 1:length(sampled_parameter_combinations_2i[1])
      push!(acts_2i,sampled_parameter_combinations_2i[1][i][1])
      push!(deacts_2i,sampled_parameter_combinations_2i[1][i][2])
      push!(degs_2i,sampled_parameter_combinations_2i[1][i][3])
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
for i in 1:length(sampled_parameter_combinations_serum[1])
      push!(acts_serum,sampled_parameter_combinations_serum[1][i][1])
      push!(deacts_serum,sampled_parameter_combinations_serum[1][i][2])
      push!(degs_serum,sampled_parameter_combinations_serum[1][i][3])
end
a=convert(Array{Float64,1},acts_serum)
b=convert(Array{Float64,1},deacts_serum)
c=convert(Array{Float64,1},degs_serum)
a2= [10^n for n in a]
b2= [10^n for n in b]
c2= [10^n for n in c]
r_deact_deg_serum =  log10.((b2)./(c2))
r_act_deg_serum =  log10.((a2)./(c2))



#paper parameters
Deg = 0.00245
Exp = 2.11
Act_serum = 0.0282
Act_2i = 0.0874
Deact_serum = 0.609
Deact_2i = 0.315
#normalise by exp
Deg = Deg/Exp
Act_serum = Act_serum/Exp
Act_2i = Act_2i/Exp
Deact_serum = Deact_serum/Exp
Deact_2i = Deact_2i/Exp
#normalise by deg
Act_norm_serum = Act_serum/Deg
Act_norm_2i = Act_2i/Deg
Deact_norm_serum  = Deact_serum/Deg
Deact_norm_2i = Deact_2i/Deg
#Transform
Act_serum_trns = [log10(Act_norm_serum)]
Act_2i_trns = [log10(Act_norm_2i)]
Deact_serum_trns = [log10(Deact_norm_serum)]
Deact_2i_trns = [log10(Deact_norm_2i)]


# fig,ax = PyPlot.subplots()
#       PyPlot.scatter(x=r_deact_deg__serum, y=r_act_deg_serum , color="cornflowerblue", label = "serum")
#       PyPlot.scatter(x=r_deact_deg_2i, y=r_act_deg_2i, color = "salmon", label = "2i")
#
#
#
# sns.jointplot


fig,ax = PyPlot.subplots()
      PyPlot.scatter([-1,2.2],[0.,0.], color="white", marker = "*")
      PyPlot.scatter([Deact_serum_trns],[Act_serum_trns], color="cornflowerblue", marker = "*")
      PyPlot.scatter([Deact_2i_trns],[Act_2i_trns], color="salmon", marker ="*")
      sns.regplot(r_deact_deg_serum,r_act_deg_serum ,color="cornflowerblue")
      sns.regplot(r_deact_deg_2i,r_act_deg_2i, color = "salmon")

      ax[:set_ylabel]("log(act/deg)")
      ax[:set_xlabel]("log(deact/deg)")
      # ax[:set_title](string(n_tester," params sampled for posterior(",n_samples,", ", tmax,") of Nanog"))
savefig(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/validate/",tau,genes[1],fileending_serum,"vs2i.png"))
