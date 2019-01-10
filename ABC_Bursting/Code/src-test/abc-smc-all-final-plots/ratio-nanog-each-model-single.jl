using PyPlot
using ABC_Bursting
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport pandas as pd

# genes = ["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]
genes = [ "ENSMUSG00000012396"]
sampled_parameter_combinations_2i = []
sampled_parameter_combinations_serum = []
m8_sampled_parameter_combinations_2i = []
m8_sampled_parameter_combinations_serum = []
m9_sampled_parameter_combinations_2i = []
m9_sampled_parameter_combinations_serum = []
n_samples = 500
n_tester = 1000
tau = 500
tmax= tau*n_samples
fileending_2i = string("_2i_0.1_fix_",n_samples,"_",tmax)
fileending_serum = string("_serum_0.1_fix_",n_samples,"_",tmax)
m8_fileending_2i = string("")
m8_fileending_serum = string("")
m9_fileending_2i = string("_20")
m9_fileending_serum = string("_20")






for name in genes
      path_to_file_2i = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/",tau,"/2i/model1_",name,fileending_2i,".txt")
      path_to_file_serum = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/",tau,"/serum/model1_",name,fileending_serum,".txt")

      m8_path_to_file_2i = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_8_0.4_5/serum_genes/2i_data/model_8_serum_genes_2i_data_", name, m8_fileending_2i,".txt")
      m8_path_to_file_serum = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_8_0.4_5/serum_genes/serum_data/model_8_serum_genes_serum_data_", name, m8_fileending_serum,".txt")
      m9_path_to_file_2i = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_9_david_02/serum_genes/2i_data/model_9_serum_genes_2i_data_", name, m9_fileending_2i,".txt")
      m9_path_to_file_serum = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_9_david_02/serum_genes/serum_data/model_9_serum_genes_serum_data_", name, m9_fileending_serum,".txt")

      data_file_2i=read_one_m1(path_to_file_2i,  3, n_samples,3)
      data_file_serum=read_one_m1(path_to_file_serum,  3, n_samples,3)

      m8_data_file_2i=read_one_m8(m8_path_to_file_2i,  3, n_samples,4)
      m8_data_file_serum=read_one_m8(m8_path_to_file_serum,  3, n_samples,4)
      m9_data_file_2i=read_one_m9(m9_path_to_file_2i,  3, n_samples,5)
      m9_data_file_serum=read_one_m9(m9_path_to_file_serum,  3, n_samples,5)

      cum_posterior_2i = Vector(n_samples)
      cum_posterior_serum = Vector(n_samples)

      m8_cum_posterior_2i = Vector(n_samples)
      m8_cum_posterior_serum = Vector(n_samples)
      m9_cum_posterior_2i = Vector(n_samples)
      m9_cum_posterior_serum = Vector(n_samples)

      cumsum!(cum_posterior_2i, data_file_2i[6])
      cumsum!(cum_posterior_serum, data_file_serum[6])

      cumsum!(m8_cum_posterior_2i, m8_data_file_2i[7])
      cumsum!(m8_cum_posterior_serum, m8_data_file_serum[7])
      cumsum!(m9_cum_posterior_2i, m9_data_file_2i[8])
      cumsum!(m9_cum_posterior_serum, m9_data_file_serum[8])

      sample_rand_2i = rand(n_tester)
      sample_rand_serum = rand(n_tester)

      m8_sample_rand_2i = rand(n_tester)
      m8_sample_rand_serum = rand(n_tester)
      m9_sample_rand_2i = rand(n_tester)
      m9_sample_rand_serum = rand(n_tester)

      winner_inds_2i = []
      winner_inds_serum = []

      m8_winner_inds_2i = []
      m8_winner_inds_serum = []
      m9_winner_inds_2i = []
      m9_winner_inds_serum = []

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

      for s in m8_sample_rand_2i
          for i in 1:length(m8_cum_posterior_2i)
                if s<m8_cum_posterior_2i[i]
                  push!(m8_winner_inds_2i, [m8_data_file_2i[1][i], m8_data_file_2i[2][i], m8_data_file_2i[3][i]])
                  break
                end
          end
      end

      for s in m8_sample_rand_serum
          for i in 1:length(m8_cum_posterior_serum)
                if s<m8_cum_posterior_serum[i]
                  push!(m8_winner_inds_serum, [m8_data_file_serum[1][i], m8_data_file_serum[2][i], m8_data_file_serum[3][i]])
                  break
                end
          end
      end
      for s in m9_sample_rand_2i
          for i in 1:length(m9_cum_posterior_2i)
                if s<m9_cum_posterior_2i[i]
                  push!(m9_winner_inds_2i, [m9_data_file_2i[1][i], m9_data_file_2i[2][i], m9_data_file_2i[3][i]])
                  break
                end
          end
      end

      for s in m9_sample_rand_serum
          for i in 1:length(m9_cum_posterior_serum)
                if s<m9_cum_posterior_serum[i]
                  push!(m9_winner_inds_serum, [m9_data_file_serum[1][i], m9_data_file_serum[2][i], m9_data_file_serum[3][i]])
                  break
                end
          end
      end


      push!(sampled_parameter_combinations_2i, winner_inds_2i)
      push!(sampled_parameter_combinations_serum, winner_inds_serum)

      push!(m8_sampled_parameter_combinations_2i, m8_winner_inds_2i)
      push!(m8_sampled_parameter_combinations_serum, m8_winner_inds_serum)
      push!(m9_sampled_parameter_combinations_2i, m9_winner_inds_2i)
      push!(m9_sampled_parameter_combinations_serum, m9_winner_inds_serum)
end



acts_2i = []
deacts_2i = []
degs_2i = []
acts_serum = []
deacts_serum = []
degs_serum = []

m8_acts_2i = []
m8_deacts_2i = []
m8_degs_2i = []
m8_acts_serum = []
m8_deacts_serum = []
m8_degs_serum = []
m9_acts_2i = []
m9_deacts_2i = []
m9_degs_2i = []
m9_acts_serum = []
m9_deacts_serum = []
m9_degs_serum = []

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

for i in 1:length(m8_sampled_parameter_combinations_2i[1])
      push!(m8_acts_2i,m8_sampled_parameter_combinations_2i[1][i][1])
      push!(m8_deacts_2i,m8_sampled_parameter_combinations_2i[1][i][2])
      push!(m8_degs_2i,m8_sampled_parameter_combinations_2i[1][i][3])
end

m8_a=convert(Array{Float64,1},m8_acts_2i)
m8_b=convert(Array{Float64,1},m8_deacts_2i)
m8_c=convert(Array{Float64,1},m8_degs_2i)
m8_a2= [10^n for n in m8_a]
m8_b2= [10^n for n in m8_b]
m8_c2= [10^n for n in m8_c]
m8_r_deact_deg_2i =  log10.((m8_b2)./(m8_c2))
m8_r_act_deg_2i =  log10.((m8_a2)./(m8_c2))

for i in 1:length(m9_sampled_parameter_combinations_2i[1])
      push!(m9_acts_2i,m9_sampled_parameter_combinations_2i[1][i][1])
      push!(m9_deacts_2i,m9_sampled_parameter_combinations_2i[1][i][2])
      push!(m9_degs_2i,m9_sampled_parameter_combinations_2i[1][i][3])
end

m9_a=convert(Array{Float64,1},m9_acts_2i)
m9_b=convert(Array{Float64,1},m9_deacts_2i)
m9_c=convert(Array{Float64,1},m9_degs_2i)
m9_a2= [10^n for n in m9_a]
m9_b2= [10^n for n in m9_b]
m9_c2= [10^n for n in m9_c]
m9_r_deact_deg_2i =  log10.((m9_b2)./(m9_c2))
m9_r_act_deg_2i =  log10.((m9_a2)./(m9_c2))




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

for i in 1:length(m8_sampled_parameter_combinations_serum[1])
      push!(m8_acts_serum, m8_sampled_parameter_combinations_serum[1][i][1])
      push!(m8_deacts_serum, m8_sampled_parameter_combinations_serum[1][i][2])
      push!(m8_degs_serum, m8_sampled_parameter_combinations_serum[1][i][3])
end
m8_a=convert(Array{Float64,1},m8_acts_serum)
m8_b=convert(Array{Float64,1},m8_deacts_serum)
m8_c=convert(Array{Float64,1},m8_degs_serum)
m8_a2= [10^n for n in m8_a]
m8_b2= [10^n for n in m8_b]
m8_c2= [10^n for n in m8_c]
m8_r_deact_deg_serum =  log10.((m8_b2)./(m8_c2))
m8_r_act_deg_serum =  log10.((m8_a2)./(m8_c2))

for i in 1:length(m9_sampled_parameter_combinations_serum[1])
      push!(m9_acts_serum, m9_sampled_parameter_combinations_serum[1][i][1])
      push!(m9_deacts_serum, m9_sampled_parameter_combinations_serum[1][i][2])
      push!(m9_degs_serum, m9_sampled_parameter_combinations_serum[1][i][3])
end
m9_a=convert(Array{Float64,1},m9_acts_serum)
m9_b=convert(Array{Float64,1},m9_deacts_serum)
m9_c=convert(Array{Float64,1},m9_degs_serum)
m9_a2= [10^n for n in m9_a]
m9_b2= [10^n for n in m9_b]
m9_c2= [10^n for n in m9_c]
m9_r_deact_deg_serum =  log10.((m9_b2)./(m9_c2))
m9_r_act_deg_serum =  log10.((m9_a2)./(m9_c2))

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

centroid_x_serum ,centroid_y_serum = (sum(r_deact_deg_serum)/length(r_deact_deg_serum), sum(r_act_deg_serum)/length(r_act_deg_serum))
centroid_x_2i ,centroid_y_2i = (sum(r_deact_deg_2i)/length(r_deact_deg_2i), sum(r_act_deg_2i)/length(r_act_deg_2i))
m8_centroid_x_serum , m8_centroid_y_serum = (sum(m8_r_deact_deg_serum)/length(m8_r_deact_deg_serum), sum(m8_r_act_deg_serum)/length(m8_r_act_deg_serum))
m8_centroid_x_2i , m8_centroid_y_2i = (sum(m8_r_deact_deg_2i)/length(m8_r_deact_deg_2i), sum(m8_r_act_deg_2i)/length(m8_r_act_deg_2i))
m9_centroid_x_serum , m9_centroid_y_serum = (sum(m9_r_deact_deg_serum)/length(m9_r_deact_deg_serum), sum(m9_r_act_deg_serum)/length(m9_r_act_deg_serum))
m9_centroid_x_2i , m9_centroid_y_2i = (sum(m9_r_deact_deg_2i)/length(m9_r_deact_deg_2i), sum(m9_r_act_deg_2i)/length(m9_r_act_deg_2i))

all_centroid_x_serum = [centroid_x_serum, m8_centroid_x_serum, m9_centroid_x_serum]
all_centroid_y_serum = [centroid_y_serum, m8_centroid_y_serum, m9_centroid_y_serum]

all_centroid_x_2i = [centroid_x_2i, m8_centroid_x_2i, m9_centroid_x_2i]
all_centroid_y_2i = [centroid_y_2i, m8_centroid_y_2i, m9_centroid_y_2i]

fig,ax = PyPlot.subplots(figsize=(5,5))
      PyPlot.plot(all_centroid_x_serum[1], all_centroid_y_serum[1], label="serum", color="cornflowerblue", markersize=30, marker = "o", linestyle="",zorder=100)
      PyPlot.plot(all_centroid_x_2i[1], all_centroid_y_2i[1], label="2i",color="salmon", marker ="o",markersize=30, linestyle="",zorder=100)
      PyPlot.plot([centroid_x_serum, centroid_x_2i], [centroid_y_serum,centroid_y_2i], "-", color="grey")
      # ax[:legend](fontsize=24, loc ="upper left")
      ax[:set_ylabel]("log(act:deg)",fontsize=24)
      ax[:yaxis][:set_label_position]("right")

      # ax[:set_xlabel]("de-act : deg",fontsize=24)
      ax[:tick_params](labelsize=18)
      ax[:set_xlim]([0.4,.9])
      ax[:set_ylim]([-.6,0.0])
      savefig(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/validate/nanog/","single_base.pdf"))

fig,ax = PyPlot.subplots(figsize=(5,5))
      PyPlot.plot(all_centroid_x_serum[2], all_centroid_y_serum[2], label="serum", color="cornflowerblue", markersize=30, marker = "o", linestyle="",zorder=100)
      PyPlot.plot(all_centroid_x_2i[2], all_centroid_y_2i[2], label="2i",color="salmon", marker ="o",markersize=30, linestyle="",zorder=100)
      PyPlot.plot([m8_centroid_x_serum, m8_centroid_x_2i], [m8_centroid_y_serum, m8_centroid_y_2i], "-", color="grey")
      # ax[:legend](fontsize=12, loc ="lower right")
      # ax[:set_ylabel]("act : deg",fontsize=24)
      # ax[:set_xlabel]("de-act : deg",fontsize=24)
      ax[:tick_params](labelsize=18)
      ax[:set_ylabel]("log(act:deg)",fontsize=24)
      ax[:yaxis][:set_label_position]("right")
      ax[:set_xlim]([-0.1,.4])
      ax[:set_ylim]([0.8,1.4])


      savefig(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/validate/nanog/","single_1.pdf"))

fig,ax = PyPlot.subplots(figsize=(5,5))
      PyPlot.plot(all_centroid_x_serum[3], all_centroid_y_serum[3], label="serum", color="cornflowerblue", markersize=30, marker = "o", linestyle="",zorder=100)
      PyPlot.plot(all_centroid_x_2i[3], all_centroid_y_2i[3], label="2i",color="salmon", marker ="o",markersize=30, linestyle="",zorder=100)
      PyPlot.plot([m9_centroid_x_serum, m9_centroid_x_2i], [m9_centroid_y_serum, m9_centroid_y_2i], "-", color="grey")
      # ax[:legend](fontsize=12, loc ="upper right")
      # ax[:set_ylabel]("act : deg",fontsize=24)
      # ax[:set_xlabel]("de-act : deg",fontsize=24)
      ax[:tick_params](labelsize=18)
      ax[:set_ylabel]("log(act:deg)",fontsize=24)
      ax[:yaxis][:set_label_position]("right")
      ax[:set_xlim]([0.0,.5])
      ax[:set_ylim]([-0.5,.1])

      savefig(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/validate/nanog/","single_2.pdf"))

fig,ax = PyPlot.subplots(figsize=(5,5))
      PyPlot.plot([Deact_2i_trns], [Act_2i_trns], color="salmon", marker ="*",  label="2i", markersize=30, linestyle="" ,zorder=100)
      PyPlot.plot([Deact_serum_trns],[Act_serum_trns], color="cornflowerblue", label="serum", markersize=30, marker = "*", linestyle="",zorder=100)
      PyPlot.plot([Deact_serum_trns, Deact_2i_trns], [Act_serum_trns, Act_2i_trns], "-", color="grey")
      # ax[:legend](fontsize=24, loc ="upper right")
      # ax[:set_ylabel]("act : deg",fontsize=24)
      # ax[:set_xlabel]("de-act : deg",fontsize=24)
      ax[:tick_params](labelsize=18)
      ax[:set_ylabel]("log(act:deg)",fontsize=24)
      ax[:yaxis][:set_label_position]("right")
      ax[:set_xlim]([2,2.5])
      ax[:set_ylim]([1,1.6])

      savefig(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/validate/nanog/","single_true.pdf"))

# print("was x 2i , y2i")
# print(Deact_2i_trns, Deact_serum_t)

#
# [0.607279, 0.137564, 0.303426]
# [-0.434196, 0.934722, -0.333503]
# was x serum , yserum
# [0.698694, 0.194088, 0.187431]
# [-0.164666, 1.26482, -0.102821]
# was x 2i , y2i
# fig,ax = PyPlot.subplots()
#
#       PyPlot.plot([Deact_serum_trns],[Act_serum_trns], color="cornflowerblue",markersize=15, marker = "*", linestyle="",zorder=100)
#       PyPlot.plot([Deact_2i_trns], [Act_2i_trns], color="salmon", marker ="*",  markersize=15, linestyle="" ,zorder=100)
#       PyPlot.plot(all_centroid_x_2i, all_centroid_y_2i, label="2i",color="salmon", marker ="o",markersize=15, linestyle="",zorder=100)
#
#        PyPlot.plot(all_centroid_x_serum, all_centroid_y_serum, label="serum", color="cornflowerblue", markersize=15, marker = "o", linestyle="",zorder=100)
#        PyPlot.plot([centroid_x_serum, centroid_x_2i], [centroid_y_serum,centroid_y_2i], "-", color="grey")
#        PyPlot.plot([m8_centroid_x_serum, m8_centroid_x_2i], [m8_centroid_y_serum, m8_centroid_y_2i], "-", color="grey")
#        PyPlot.plot([m9_centroid_x_serum, m9_centroid_x_2i], [m9_centroid_y_serum, m9_centroid_y_2i], "-", color="grey")
#        PyPlot.plot([Deact_serum_trns, Deact_2i_trns], [Act_serum_trns, Act_2i_trns], "-", color="grey")
#
#        ax[:annotate]("Base",xy=(centroid_x_2i, centroid_y_2i), xytext = (centroid_x_2i+0.1, centroid_y_2i+0.1),fontsize = 12,zorder=100)
#        ax[:annotate]("Feedback 1",xy=(m8_centroid_x_2i, m8_centroid_y_2i), xytext= (m8_centroid_x_2i+0.1, m8_centroid_y_2i+0.1), fontsize = 12,zorder=100)
#        ax[:annotate]("Feedback 2",xy=(m9_centroid_x_2i, m9_centroid_y_2i), xytext= (m9_centroid_x_2i+0.1, m9_centroid_y_2i+0.1), fontsize = 12,zorder=100)
#       ax[:legend](fontsize=12, loc ="lower right")
#       ax[:set_ylabel]("log(act/deg)",fontsize=12)
#       ax[:set_xlabel]("log(deact/deg)",fontsize=12)
#       ax[:tick_params](labelsize=12)
#       savefig(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/validate/nanog/","all.pdf"))
