using PyCall, PyPlot; @pyimport pandas as pan
@pyimport seaborn as sns
using ABC_Bursting
using KernelDensity

name="ENSMUSG00000003032"
a3_2i,b3_2i,c3_2i,d3_2i,di3_2i,w3_2i = 0,0,0,0,0,0
fileending_2i="validate_model_9_david_100_01_250_3_0c2_no_k"
nsamples = 100
path_to_file_2i = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9/",fileending_2i,".txt")
path_to_plot = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/validation-simulation/model-9/no_deact_k/",name,fileending_2i,"_joint_validate.pdf")
a3_2i,b3_2i,c3_2i,d3_2i,di3_2i,w3_2i=read_one_m1(path_to_file_2i,  3, nsamples, 3)
d_2i = hcat(convert(Array{Float64,1}, a3_2i),
              convert(Array{Float64,1}, b3_2i),
              convert(Array{Float64,1}, c3_2i))
a_2i = pan.DataFrame(data= Dict( :c_act=>d_2i[:,1],:c_deact=>d_2i[:,2],:c_deg=>d_2i[:,3]))

### calculate joint approx #####
a3_2i_float =convert(Array{Float64}, a3_2i)
b3_2i_float =convert(Array{Float64}, b3_2i)
c3_2i_float =convert(Array{Float64}, c3_2i)
a3_2i_float_k = KernelDensity.kde(a3_2i_float)
b3_2i_float_k = KernelDensity.kde(b3_2i_float)
c3_2i_float_k = KernelDensity.kde(c3_2i_float)
pa = []
pb = []
pc = []
for i in 1:length(a3_2i)
              p_a3_2i = pdf(a3_2i_float_k, a3_2i[i])
              p_b3_2i = pdf(b3_2i_float_k, b3_2i[i])
              p_c3_2i = pdf(c3_2i_float_k, c3_2i[i])
              push!(pa, p_a3_2i)
              push!(pb, p_b3_2i)
              push!(pc, p_c3_2i)
end
all_a=sum(pa)
all_b=sum(pb)
all_c=sum(pc)
pa_norm = pa./all_a
pb_norm = pb./all_b
pc_norm = pc./all_c
p=[]
for i in 1:length(a3_2i)
              pi = pa_norm[i]*pb_norm[i]* pc_norm[i]
              push!(p,pi)
end
maxe=indmax(p)
###
PyPlot.ion()
fig, ax = PyPlot.subplots(figsize = (15,15), ncols=3, nrows=3)
      PyPlot.ioff()
      cmap = sns.cubehelix_palette(as_cmap=true, dark=0, light=1, reverse=true)
      newTrue=log10.([0.001, 0.0001, 0.04])
      sns.distplot(a3_2i,bins=15, kde= true, hist=true,rug= true ,color = "grey",ax=ax[1,1])
      sns.distplot(b3_2i,bins=15, kde= true, hist=true,rug= true ,   color = "grey",ax=ax[2,2])
      sns.distplot(c3_2i,bins=15, kde= true, hist=true,rug= true ,  color = "grey", ax=ax[3,3])
      ax[2,1][:set_xlim]([min(a3_2i...)-0.31,max(a3_2i...)+0.3])
      ax[2,1][:set_ylim]([min(b3_2i...),max(b3_2i...)])
      ax[3,1][:set_xlim]([min(a3_2i...)-0.31,max(a3_2i...)+0.3])
      ax[3,1][:set_ylim]([min(c3_2i...),max(c3_2i...)+0.3])
      ax[3,2][:set_xlim]([min(b3_2i...)-0.2,max(b3_2i...)+0.2])
      ax[3,2][:set_ylim]([min(c3_2i...),max(c3_2i...)+0.3])


      sns.jointplot(x="c_act", y="c_deact", data=a_2i, kind="kde",color="gray", ax=ax[2,1])
      sns.jointplot(x="c_act", y="c_deg", data=a_2i, kind="kde", color="grey",ax=ax[3,1])
      sns.jointplot(x="c_deact", y="c_deg", data=a_2i, kind="kde", color="grey",ax=ax[3,2])

      ax[1,1][:set_ylabel]("Activation rate",fontsize=19)
      ax[1,1][:tick_params](labelsize=12)
      # ax[2,2][:set_xlabel]("Deactivation rate",fontsize=12)
      ax[2,2][:tick_params](labelsize=12)
      ax[3,3][:set_xlabel]("Base line",fontsize=19)
      ax[3,3][:tick_params](labelsize=12)

      ax[2,1][:set_ylabel]("Degradation rate",fontsize=19)
      # ax[2,1][:set_xlabel]("Activation rate",fontsize=12)
      ax[2,1][:tick_params](labelsize=12)

      ax[3,1][:set_ylabel]("Base line",fontsize=19)
      ax[3,1][:set_xlabel]("Activation rate",fontsize=19)
      ax[3,1][:tick_params](labelsize=12)

      # ax[3,2][:set_ylabel]("Deactivation rate",fontsize=12)
      ax[3,2][:set_xlabel]("Degradation rate",fontsize=19)
      ax[3,2][:tick_params](labelsize=12)



      ax[2,1][:scatter]([truth[1]], [truth[2]],color ="gold", marker = "x",linewidth=3,  s =600, zorder =120)
      ax[3,1][:scatter]([truth[1]], [truth[3]],color ="gold", marker = "x",linewidth=3,  s =600, zorder =120)
      ax[3,2][:scatter]([truth[2]], [truth[3]],color ="gold", marker = "x",linewidth=3,  s =600, zorder =120)

      ax[1,1][:axvline]([truth[1]],color ="gold", label = "true value",linewidth=3)
      ax[2,2][:axvline]([truth[2]],color ="gold",linewidth=3)
      ax[3,3][:axvline]([truth[3]],color ="gold",linewidth=3)

      ax[1,1][:axvline](a3_2i[maxe],color ="purple",linewidth=3, label = "max. joint probability")
      ax[2,2][:axvline](b3_2i[maxe],color ="purple",linewidth=3)
      ax[3,3][:axvline](c3_2i[maxe],color ="purple",linewidth=3)
      ax[1,1][:legend](fontsize = 18, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=1.8)
      # ax[1,1][:axvline](a3_2i[indmin(di3_2i)],color ="green")
      # ax[2,2][:axvline](b3_2i[indmin(di3_2i)],color ="green")
      # ax[3,3][:axvline](c3_2i[indmin(di3_2i)],color ="green")
      #
      # ax[1,1][:axvline](mean(a3_2i),color ="pink")
      # ax[2,2][:axvline](mean(b3_2i),color ="pink")
      # ax[3,3][:axvline](mean(c3_2i),color ="pink")

      ax[1,2][:axis]("off")
      ax[2,3][:axis]("off")
      ax[1,3][:axis]("off")


      # ax[1,3][:annotate]("gold:   true value",xy=(.1,.1), xytext= (.2,.4),color="gold", fontsize = 19,zorder=100)
      # ax[1,3][:annotate]("purple: maximal joint probability",xy=(.1,.1), xytext= (.2,.6),color="purple", fontsize = 19,zorder=100)
      # ax[1,3][:annotate]("green:min distance",xy=(.1,.1), xytext= (.2,.8),color="green", fontsize = 19,zorder=100)
      # ax[1,3][:annotate]("pink: mean marginal",xy=(.1,.1), xytext= (.2,.2),color="pink", fontsize = 19,zorder=100)
      # ax[1,3][:tick_params](color="white",labelsize=1)
      # ax[1,3][:spines]["top"][:set_color]("none")
      # ax[1,3][:spines]["right"][:set_color]("none")
      # ax[1,3][:spines]["left"][:set_color]("none")
      # ax[1,3][:spines]["bottom"][:set_color]("none")


PyPlot.savefig(path_to_plot)

# spines['bottom'].set_color('#dddddd')
# ax.yaxis.label.set_color('red')
