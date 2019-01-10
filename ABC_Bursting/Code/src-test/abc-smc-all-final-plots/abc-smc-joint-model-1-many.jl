using PyCall, PyPlot; @pyimport pandas as pd
using PyCall, PyPlot; @pyimport seaborn as sns
using ABC_Bursting
using StatPlots
#names=["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]

names=["ENSMUSG00000012396"]
#,"ENSMUSG00000027715"]
#, "ENSMUSG00000029177"]
#names=["ENSMUSG00000006398"]
a3_2i,b3_2i,c3_2i,d3_2i,di3_2i,w3_2i = 0,0,0,0,0,0
a3_serum,b3_serum,c3_serum,d3_serum,di3_serum,w3_serum = 0,0,0,0,0,0


fileending_2i="_2i_0.1_fix_500_250000"
fileending_serum="_serum_0.1_fix_500_250000"
p_1 = [-6,0.]
# p_2 = [-6,0]
# p_3 = [-6,0]

nsamples = 500
for name in names

  path_to_file_2i = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/500/2i/model1_",name,fileending_2i,".txt")
  path_to_file_serum = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/500/serum/model1_",name,fileending_serum,".txt")
  path_to_plot = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/model-1-nanog/500/serum/",name,fileending_2i,fileending_serum,"_joint.png")
  a3_2i,b3_2i,c3_2i,d3_2i,di3_2i,w3_2i=read_one_m1(path_to_file_2i,  3, nsamples, 3)
  a3_serum,b3_serum,c3_serum,d3_serum,di3_serum,w3_serum=read_one_m1(path_to_file_serum,  3, nsamples, 3)

  d_2i = hcat(convert(Array{Float64,1}, a3_2i),
              convert(Array{Float64,1}, b3_2i),
              convert(Array{Float64,1}, c3_2i))

  d_serum = hcat(convert(Array{Float64,1}, a3_serum),
              convert(Array{Float64,1}, b3_serum),
              convert(Array{Float64,1}, c3_serum))
  PyPlot.ion()

  a_2i = pd.DataFrame(data= Dict( :c_act=>d_2i[:,1],:c_deact=>d_2i[:,2],:c_deg=>d_2i[:,3]))
  a_serum = pd.DataFrame(data= Dict( :c_act=>d_serum[:,1],:c_deact=>d_serum[:,2],:c_deg=>d_serum[:,3]))

      fig, ax = PyPlot.subplots(figsize=(15,15),ncols=3, nrows=3)
      left   =  0.125  # the left side of the subplots of the figure
      right  =  0.9    # the right side of the subplots of the figure
      bottom =  0.1    # the bottom of the subplots of the figure
      top    =  0.9    # the top of the subplots of the figure
      wspace =  .5     # the amount of width reserved for blank space between subplots
      hspace =  .5    # the amount of height reserved for white space between subplots
      PyPlot.subplots_adjust(
          left    =  left,
          bottom  =  bottom,
          right   =  right,
          top     =  top,
          wspace  =  wspace,
          hspace  =  hspace
      )
      PyPlot.ioff()
      cmap = sns.cubehelix_palette(as_cmap=true, dark=0, light=1, reverse=true)



      sns.distplot(a3_serum,bins=15, kde= true, hist=true,rug= false ,color = "cornflowerblue", label = "lisi",axlabel = "c act",ax=ax[1,1])
      sns.distplot(b3_serum,bins=15, kde= true, hist=true,rug= false ,   color = "cornflowerblue", label = "lisi",axlabel = "c de-act",ax=ax[2,2])
      sns.distplot(c3_serum,bins=15, kde= true, hist=true,rug= false ,  color = "cornflowerblue", label = "lisi",axlabel = "c deg",ax=ax[3,3])
      sns.distplot(a3_2i,bins=15, kde= true, hist=true,rug= false ,   color = "salmon", label = "lisi",axlabel = "Activation rate cact",ax=ax[1,1])
     sns.distplot(b3_2i,bins=15, kde= true, hist=true,rug= false ,  color = "salmon", label = "lisi",axlabel = "Deactivation rate cde-act",ax=ax[2,2])
     sns.distplot(c3_2i,bins=15, kde= true, hist=true,rug= false ,   color = "salmon", label = "lisi",axlabel = "degradation rate cdeg",ax=ax[3,3])
      sns.jointplot(x="c_act", y="c_deact", data=a_2i, kind="kde",color="lightpink", ax=ax[2,1])
      sns.jointplot(x="c_act", y="c_deg", data=a_2i, kind="kde", color="lightpink",ax=ax[3,1])
      sns.jointplot(x="c_deact", y="c_act", data=a_2i, kind="kde", color="lightpink",ax=ax[3,2])

      sns.jointplot(y="c_act", x="c_deact", data=a_serum, kind="kde", color="cornflowerblue", ax=ax[1,2])
      sns.jointplot(y="c_act", x="c_deg", data=a_serum, kind="kde", color="cornflowerblue",ax=ax[1,3])
      sns.jointplot(y="c_deact", x="c_act", data=a_serum, kind="kde", color="cornflowerblue",ax=ax[2,3])

      PyPlot.show(fig)
      PyPlot.savefig(path_to_plot)
end
