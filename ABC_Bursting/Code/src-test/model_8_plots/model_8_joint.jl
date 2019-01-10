using PyCall, PyPlot; @pyimport pandas as pd
using PyCall, PyPlot; @pyimport seaborn as sns
using ABC_Bursting
using StatPlots
# genes=["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867"] #5 first
# , "ENSMUSG00000041431"]
# genes = ["ENSMUSG00000003032", "ENSMUSG00000012396", "ENSMUSG00000018604", "ENSMUSG00000020717", "ENSMUSG00000021255", "ENSMUSG00000021835", "ENSMUSG00000022528", "ENSMUSG00000025056", "ENSMUSG00000038793", "ENSMUSG00000048402", "ENSMUSG00000051176"]
# genes = ["ENSMUSG00000003032", "ENSMUSG00000012396", "ENSMUSG00000018604", "ENSMUSG00000020717", "ENSMUSG00000021255"]  #5 fisrt serum
# , "ENSMUSG00000021835", "ENSMUSG00000022528", "ENSMUSG00000025056", "ENSMUSG00000038793", "ENSMUSG00000048402", "ENSMUSG00000051176"]

#,"ENSMUSG00000027715"]
#, "ENSMUSG00000029177"]
genes=["ENSMUSG00000012396"]
a3_2i,b3_2i,c3_2i,d3_2i,di3_2i,w3_2i = 0,0,0,0,0,0
a3_serum,b3_serum,c3_serum,d3_serum,di3_serum,w3_serum = 0,0,0,0,0,0


fileending_2i="model_8_serum_genes_2i_data_"
fileending_serum="model_8_serum_genes_serum_data_"
p_1 = [-6,0.]
# p_2 = [-6,0]
# p_3 = [-6,0]

nsamples = 500
for name in genes
  path_to_file_2i = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_8_02/serum_genes/2i_data/",fileending_2i,name,".txt")
  path_to_file_serum = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_8_02/serum_genes/serum_data/",fileending_serum,name,".txt")
  path_to_plot = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/model-8-joint-02/",name,fileending_serum,"joint1.pdf")
  print("path_to_file_2i", path_to_file_2i,"\n")
  print("path_to_file_serum", path_to_file_serum,"\n")
  a3_2i,b3_2i,c3_2i,d3_2i,e3_2i,di3_2i,w3_2i=read_one_m8(path_to_file_2i,  3, nsamples, 4)
  a3_serum,b3_serum,c3_serum,d3_serum,e3_serum,di3_serum,w3_serum=read_one_m8(path_to_file_serum,  3, nsamples, 4)

  d_2i = hcat(convert(Array{Float64,1}, a3_2i),
              convert(Array{Float64,1}, b3_2i),
              convert(Array{Float64,1}, c3_2i),
              convert(Array{Float64,1}, d3_2i))

  d_serum = hcat(convert(Array{Float64,1}, a3_serum),
              convert(Array{Float64,1}, b3_serum),
              convert(Array{Float64,1}, c3_serum),
              convert(Array{Float64,1}, d3_serum))
  PyPlot.ion()

  a_2i = pd.DataFrame(data= Dict( :c_act=>d_2i[:,1],:c_deact=>d_2i[:,2],:c_deg=>d_2i[:,3],:fb=>d_2i[:,4]))
  a_serum = pd.DataFrame(data= Dict( :c_act=>d_serum[:,1],:c_deact=>d_serum[:,2],:c_deg=>d_serum[:,3],:fb=>d_serum[:,4]))

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
      cmap = sns.cubehelix_palette(as_cmap=true, dark=0, light=1, reverse=true)



      sns.distplot(a3_serum,bins=15, kde= true, hist=true,rug= false ,color = "cornflowerblue", label = "lisi",ax=ax[1,1])
      sns.distplot(b3_serum,bins=15, kde= true, hist=true,rug= false ,   color = "cornflowerblue", label = "lisi",ax=ax[2,2])
      sns.distplot(c3_serum,bins=15, kde= true, hist=true,rug= false ,  color = "cornflowerblue", label = "lisi",ax=ax[3,3])
      sns.distplot(d3_serum,bins=15, kde= true, hist=true,rug= false ,  color = "cornflowerblue", label = "lisi",ax=ax[4,4])

      sns.distplot(a3_2i,bins=15, kde= true, hist=true,rug= false ,   color = "salmon", label = "lisi",ax=ax[1,1])
     sns.distplot(b3_2i,bins=15, kde= true, hist=true,rug= false ,  color = "salmon", label = "lisi",ax=ax[2,2])
     sns.distplot(c3_2i,bins=15, kde= true, hist=true,rug= false ,   color = "salmon", label = "lisi",ax=ax[3,3])
     sns.distplot(d3_2i,bins=15, kde= true, hist=true,rug= false ,   color = "salmon", label = "lisi",ax=ax[4,4])

     sns.jointplot(x="c_act", y="c_deact", data=a_2i, kind="kde",color="lightpink", ax=ax[2,1])
     sns.jointplot(x="c_act", y="c_deg", data=a_2i, kind="kde", color="lightpink",ax=ax[3,1])
     sns.jointplot(x="c_act", y="fb", data=a_2i, kind="kde", color="lightpink",ax=ax[4,1])
     sns.jointplot(x="c_deact", y="c_deg", data=a_2i, kind="kde", color="lightpink",ax=ax[3,2])
     sns.jointplot(x="c_deact", y="fb", data=a_2i, kind="kde", color="lightpink",ax=ax[4,2])
     sns.jointplot(x="c_deg", y="fb", data=a_2i, kind="kde", color="lightpink",ax=ax[4,3])


     sns.jointplot(x="c_act", y="c_deact", data=a_serum, kind="kde", color="cornflowerblue", ax=ax[1,2])
     sns.jointplot(x="c_act", y="c_deg", data=a_serum, kind="kde", color="cornflowerblue",ax=ax[1,3])
     sns.jointplot(x="c_act", y="fb", data=a_serum, kind="kde", color="cornflowerblue",ax=ax[1,4])
     sns.jointplot(x="c_deact", y="c_deg", data=a_serum, kind="kde", color="cornflowerblue",ax=ax[2,3])
     sns.jointplot(x="c_deact", y="fb", data=a_serum, kind="kde", color="cornflowerblue",ax=ax[2,4])
     sns.jointplot(x="c_deg", y="fb", data=a_serum, kind="kde", color="cornflowerblue",ax=ax[3,4])
###

     ax[1,1][:set_ylabel]("Activation rate",fontsize=16)
     ax[2,1][:set_ylabel]("Deactivation rate",fontsize=16)
     ax[3,1][:set_ylabel]("Degradation rate",fontsize=16)
     ax[4,1][:set_ylabel]("Feedback k",fontsize=16)
     ax[4,1][:set_xlabel]("Activation rate",fontsize=16)
     ax[4,2][:set_xlabel]("Deactivation rate",fontsize=16)
     ax[4,3][:set_xlabel]("Degradation rate",fontsize=16)
     ax[4,4][:set_xlabel]("Feedback k",fontsize=16)

    ax[1,1][:tick_params](labelsize=12)
    ax[2,2][:tick_params](labelsize=12)
    ax[3,3][:tick_params](labelsize=12)
    ax[4,4][:tick_params](labelsize=12)
    ax[2,1][:tick_params](labelsize=12)
    ax[3,1][:tick_params](labelsize=12)
    ax[4,1][:tick_params](labelsize=12)
    ax[3,2][:tick_params](labelsize=12)
    ax[4,2][:tick_params](labelsize=12)
    ax[4,3][:tick_params](labelsize=12)

    ax[1,2][:tick_params](labelsize=12)
    ax[1,3][:tick_params](labelsize=12)
    ax[1,4][:tick_params](labelsize=12)
    ax[2,3][:tick_params](labelsize=12)
    ax[2,4][:tick_params](labelsize=12)
    ax[3,4][:tick_params](labelsize=12)
    ##
      PyPlot.show(fig)
      PyPlot.savefig(path_to_plot)
end
      #

      # sns.jointplot(x="c_act", y="c_deact", data=a_2i, kind="kde",color="lightpink", ax=ax[2,1])
      # sns.jointplot(x="c_act", y="c_deg", data=a_2i, kind="kde", color="lightpink",ax=ax[3,1])
      # sns.jointplot(x="c_act", y="fb", data=a_2i, kind="kde", color="lightpink",ax=ax[4,1])
      # sns.jointplot(x="c_deact", y="c_deg", data=a_2i, kind="kde", color="lightpink",ax=ax[3,2])
      # sns.jointplot(x="c_deact", y="fb", data=a_2i, kind="kde", color="lightpink",ax=ax[4,2])
      # sns.jointplot(x="c_deg", y="fb", data=a_2i, kind="kde", color="lightpink",ax=ax[4,3])
      #
      #
      # sns.jointplot(y="c_act", x="c_deact", data=a_serum, kind="kde", color="cornflowerblue", ax=ax[1,2])
      # sns.jointplot(y="c_act", x="c_deg", data=a_serum, kind="kde", color="cornflowerblue",ax=ax[1,3])
      # sns.jointplot(y="c_act", x="fb", data=a_serum, kind="kde", color="cornflowerblue",ax=ax[1,4])
      # sns.jointplot(y="c_deact", x="c_deg", data=a_serum, kind="kde", color="cornflowerblue",ax=ax[2,3])
      # sns.jointplot(y="c_deact", x="fb", data=a_serum, kind="kde", color="cornflowerblue",ax=ax[2,4])
      # sns.jointplot(y="c_deg", x="fb", data=a_serum, kind="kde", color="cornflowerblue",ax=ax[3,4])
