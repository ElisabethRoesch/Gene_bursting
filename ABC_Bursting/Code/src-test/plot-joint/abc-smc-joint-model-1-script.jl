using PyCall, PyPlot; @pyimport pandas as pd
using PyCall, PyPlot; @pyimport seaborn as sns

using ABC_Bursting
using StatPlots
#names=["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]

names=["ENSMUSG00000012396"]
#,"ENSMUSG00000027715"]
#, "ENSMUSG00000029177"]
#names=["ENSMUSG00000006398"]
a3,b3,c3,d3,di3,w3 = 0,0,0,0,0,0
fileending="_serum_0.1_fix_500_250000"
nsamples = 500
for name in names
    #model1_ENSMUSG00000006398_2i_0.1_fix_100
  pathtofile_lisi = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/500/serum/model1_",name,fileending,".txt")
  # pathtoplotA = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/testing-priors/2_6_40_80/",name,fileending,"_plotA.png")
  pathtoplotB = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/model-1-nanog/500/serum/",name,fileending,"_joint.png")
  a3,b3,c3,d3,di3,w3=read_one_m1(pathtofile_lisi,  3, nsamples, 3)
  d=hcat(convert(Array{Float64,1}, a3),
      convert(Array{Float64,1}, b3),
      convert(Array{Float64,1}, c3))

      PyPlot.ion()

  a = pd.DataFrame(data= Dict( :c_act=>d[:,1],:c_deact=>d[:,2],:c_deg=>d[:,3]))
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
      sns.distplot(a3,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c act",ax=ax[1,1])
      sns.distplot(b3,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c de-act",ax=ax[2,2])
      sns.distplot(c3,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c deg",ax=ax[3,3])
      # sns.distplot(d3,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "scaling factor",ax=ax[4,4])
      sns.jointplot(x="c_act", y="c_deact", data=a, kind="kde", color="m", ax=ax[2,1])
      # sns.kdeplot(a[:c_act], a[:c_deg],  cmap=cmap, n_levels=60, shade=true, ax=ax[2,1])
      sns.jointplot(x="c_act", y="c_deg", data=a, kind="kde", color="m",ax=ax[3,1])
      sns.jointplot(x="c_deact", y="c_act", data=a, kind="kde", color="m",ax=ax[3,2])
      # sns.jointplot(x="c_act", y="scale", data=a, kind="kde", color="m",ax=ax[4,1])
      # sns.jointplot(x="c_deact", y="scale", data=a, kind="kde", color="m",ax=ax[4,2])
      # sns.jointplot(x="c_act", y="scale", data=a, kind="kde", color="m",ax=ax[4,3])
      sns.distplot(a3,bins=20, kde= true, hist=true,rug= true ,  color = "white", label = "lisi",axlabel = "c act",ax=ax[1,2])
      PyPlot.show(fig)
      PyPlot.savefig(pathtoplotB)
end

PyPlot.show(fig)
PyPlot.savefig(pathtoplotB)


PyPlot.ion()
a = pd.DataFrame(data= Dict( :c_act=>d[:,1],:c_deact=>d[:,2],:c_deg=>d[:,3],:scale=>d[:,4]))

a1=sns.jointplot(x="c_act", y="scale", data=a, kind="kde", color="m")
b= sns.jointplot(x="c_deact", y="scale", data=a, kind="kde", color="m")
show(a1)
savefig("me.png")

;
d=hcat(convert(Array{Float64,1}, a3),
    convert(Array{Float64,1}, b3),
    convert(Array{Float64,1}, c3),
    convert(Array{Float64,1}, d3))
StatPlots.corrplot(d)
a = pd.DataFrame(data= Dict( :c_act=>d[:,1],:c_deact=>d[:,2],:c_deg=>d[:,3],:scale=>d[:,4]))
PyPlot.ion()

g = sns.jointplot(x="c_act", y="scale", data=a, kind="kde", color="m")
g = sns.jointplot(x="c_deact", y="scale", data=a, kind="kde", color="m")
g = sns.jointplot(x="c_deact", y="scale", data=a, kind="kde", color="m")
# g = sns.jointplot(x="c_act", y="scale", data=a, kind="kde", color="m")
g = sns.PairGrid(a)
sns.jointplot(x="c_act", y="scale", data=a, kind="kde", color="m")
g.plot( sns.jointplot(x="c_act", y="scale", data=a))

PyPlot.ion()
fig, ax = PyPlot.subplots(figsize=(15,15),ncols=4, nrows=4)
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
sns.distplot(a3,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c act",ax=ax[1,1])
sns.distplot(b3,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c de-act",ax=ax[2,2])
sns.distplot(c3,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c deg",ax=ax[3,3])
sns.distplot(d3,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "scaling factor",ax=ax[4,4])
sns.jointplot(x="c_act", y="c_deact", data=a, kind="kde", color="m", ax=ax[2,1])
sns.jointplot(x="c_act", y="c_deg", data=a, kind="kde", color="m",ax=ax[3,1])
sns.jointplot(x="c_deact", y="c_act", data=a, kind="kde", color="m",ax=ax[3,2])
sns.jointplot(x="c_act", y="scale", data=a, kind="kde", color="m",ax=ax[4,1])
sns.jointplot(x="c_deact", y="scale", data=a, kind="kde", color="m",ax=ax[4,2])
sns.jointplot(x="c_act", y="scale", data=a, kind="kde", color="m",ax=ax[4,3])
PyPlot.suptitle("Original vs Normalized vs Standardized")
