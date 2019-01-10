using PyCall, PyPlot; @pyimport seaborn as sns
using ABC_Bursting

# names=["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]

names=["ENSMUSG00000041431"]
# , "ENSMUSG00000029177"]
#model1_ENSMUSG00000041431_0.1_500_250000_2i_genes_2i_data.txt
 a,b,c,a2,b2,c2,a3,b3,c3 =[],[],[],[],[],[],[],[],[]
# /project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/explore_distr/model_1/all_high_var_in_2i_genes_in_2
for name in names
  pathtofile_lisi = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/explore_distr/model_1/all_high_var_in_2i_genes_in_2i/model1_",name,"_0.1_500_250000_2i_genes_2i_data.txt")

  # pathtofile_lisi = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/old/model_1_2i_medium/model1_",name,"_2i_0.1_100.txt")
  pathtoplot = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/old/model-1-medium/old",name)
   a,b,c,a2,b2,c2,a3,b3,c3 = plot_one_abacus_file_m1(pathtofile_lisi, pathtoplot, 3, 100, 4)
end



old_prior = [-3,3]
new_prior = [-6,0]
PyPlot.ion()
fig, ax = PyPlot.subplots(figsize=(15,5),ncols=3, nrows=1)
   left   =  0.125  # the left side of the subplots of the figure
   right  =  0.9    # the right side of the subplots of the figure
   bottom =  0.1    # the bottom of the subplots of the figure
   top    =  0.9    # the top of the subplots of the figure
   wspace =  .5     # the amount of width reserved for blank space between subplots
   hspace =  1.1    # the amount of height reserved for white space between subplots
   PyPlot.subplots_adjust(
       left    =  left,
       bottom  =  bottom,
       right   =  right,
       top     =  top,
       wspace  =  wspace,
       hspace  =  hspace
   )
   y_title_margin = 1.2
   PyPlot.suptitle("Original vs Normalized vs Standardized", y = 1.09, fontsize=20)
   sns.distplot(b,bins=20, kde= true, hist=true,rug= true ,  color = "green",  label = "epsilon = 0.4",ax=ax[1,1])
   sns.distplot(b2,bins=20, kde= true, hist=true,rug= true ,  color = "green",  label = "epsilon = 0.2",ax=ax[2,1])
   sns.distplot(b3,bins=20, kde= true, hist=true,rug= true ,  color = "green",label = "epsilon = 0.1",ax=ax[3,1])

   ax[1,1][:axvline]([new_prior[1]],color ="orange",linewidth=3,label = "prior")
   ax[2,1][:axvline]([new_prior[1]],color ="orange",linewidth=3,label = "prior")
   ax[3,1][:axvline]([new_prior[1]],color ="orange",linewidth=3,label = "prior")
   ax[1,1][:axvline]([new_prior[2]],color ="orange",linewidth=3)
   ax[2,1][:axvline]([new_prior[2]],color ="orange",linewidth=3)
   ax[3,1][:axvline]([new_prior[2]],color ="orange",linewidth=3)
   ax[1,1][:tick_params](labelsize=14)
   ax[2,1][:tick_params](labelsize=14)
   ax[3,1][:tick_params](labelsize=14)
   ax[1,1][:set_ylabel]("New prior",fontsize=19)
   ax[1,1][:set_xlabel]("Deactivation rate",fontsize=12)
   ax[2,1][:set_xlabel]("Deactivation rate",fontsize=12)
   ax[3,1][:set_xlabel]("Deactivation rate",fontsize=12)
   ax[1,1][:set_xlim]([new_prior[1]-1,old_prior[2]+1])
   ax[2,1][:set_xlim]([new_prior[1]-1,old_prior[2]+1])
   ax[3,1][:set_xlim]([new_prior[1]-1,old_prior[2]+1])
   ax[1,1][:legend](fontsize=12)
   ax[2,1][:legend](fontsize=12)
   ax[3,1][:legend](fontsize=12)

 PyPlot.savefig("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/old/model-1-medium/new.pdf")
