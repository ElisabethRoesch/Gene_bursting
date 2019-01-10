using PyPlot
using ABC_Bursting
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport pandas as pd

# genes = ["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]
genes = [ "ENSMUSG00000048402","ENSMUSG00000025056", "ENSMUSG00000048402", "ENSMUSG00000051176"]
colors = ["darkblue", "cornflowerblue", "orangered", "darkred"]

real_names = ["Gli2 in serum","Nr0b1 in serum",  "Gli2 in 2i" , "Zfp42 in 2i" ]
# markers = ["+" , "*" , "^","+" , "*" , "^"]
sampled_parameter_combinations_2i = []
sampled_parameter_combinations_serum = []
n_samples = 500
n_tester = 100

all_r_deact_deg__serum = []
all_r_act_deg_serumm = []

endfile = ["serum","serum","2i","2i"]
counter = 0
for name in genes
      counter=counter+1

      path_to_file_serum =  string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/explore_distr/model_1/",endfile[counter],"/model1_",name, "_",endfile[counter],"_0.1_fix_500_250000.txt")
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
      r_deact_deg__serum =  log10.((b2)./(c2))
      r_act_deg_serum =  log10.((a2)./(c2))
      push!(all_r_deact_deg__serum, r_deact_deg__serum)
      push!(all_r_act_deg_serumm, r_act_deg_serum)
end


fig,ax = PyPlot.subplots()
      ax[:axvline]([0],color ="black",alpha=0.3, linestyle ="--")
      ax[:axhline]([0],color ="black",alpha=0.3, linestyle ="--")
      sns.regplot([-2.5,2.5],[0.0,-2.0],color = "white",fit_reg=false)
      for name in 1:length(genes)
            sns.regplot(x=all_r_deact_deg__serum[length(genes)-name+1],y=all_r_act_deg_serumm[length(genes)-name+1],label=string(real_names[length(genes)-name+1]),color = colors[length(genes)-name+1], fit_reg=false,scatter = true)
      end
      ax[:legend](fontsize=11)
      ax[:set_ylabel]("log(act/deg)",fontsize=16)
      ax[:set_xlabel]("log(deact/deg)",fontsize=16)
      ax[:tick_params](labelsize=12)
      ax[:annotate]("I.",xy=(.1,.1), xytext= (2.4,.4), color ="grey",fontsize = 19,zorder=100)
      ax[:annotate]("II.",xy=(.1,.1), xytext= (-2.5,0.4), color ="grey",fontsize = 19,zorder=100)
      ax[:annotate]("III.",xy=(.1,.1), xytext= (-2.6,-.2), color ="grey",fontsize = 19,zorder=100)
      ax[:annotate]("VI.",xy=(.1,.1), xytext= (2.3,-.2),color ="grey", fontsize = 19,zorder=100)



      # ax[:set_title](string(" params sampled for ", genes,))
      # ax[:set_ylim]([-2.5,1])
      # ax[:set_xlim]([-2.5,2.5])
PyPlot.savefig(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/validate/no_line_four_log_zerolines.pdf"))
