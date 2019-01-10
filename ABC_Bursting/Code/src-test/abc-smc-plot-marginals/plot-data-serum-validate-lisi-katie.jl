using PyCall, PyPlot; @pyimport seaborn as sns

pathtofile_lisi = ""
pathtofile_katie = ""
pathtoplot= ""
if length(ARGS) == 3
    print("you have passed 3 arguments. i am thinking they are lisis estimation, katies estimation and path to plot")
    pathtofile_lisi = ARGS[1]
    pathtofile_katie = ARGS[2]
    pathtoplot= ARGS[3]
elseif length(ARGS) == 2
    print("you have passed 2 arguments. i am thinking they are gene name and path to plot")
    name = ARGS[1]
    pathtofile_lisi = string(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/lisi_high_var_serum/model1_",name),"_serum_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt")
    pathtofile_katie = string(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/ABCout_RNAseq_transformed_keygenes_2_serum/model1_",name),"_serum_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt")
    pathtoplot= ARGS[2]
elseif length(ARGS) == 1
    print("you have passed 1 argument. i am thinking it is gene name ")
    name = ARGS[1]
    pathtofile_lisi = string(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/lisi_high_var_serum/model1_",name),"_serum_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt")
    pathtofile_katie = string(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/ABCout_RNAseq_transformed_keygenes_2_serum/model1_",name),"_serum_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt")
    pathtoplot = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/",name)
end



# pathtofile_lisi ="/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/lisi_high_var_serum/model1_ENSMUSG00000003032_serum_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt"
# pathtofile_katie ="/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/ABCout_RNAseq_transformed_keygenes_2_serum/model1_ENSMUSG00000003032_serum_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt"
# pathtoplot= "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/a.pdf"

# pathtofile ="/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/ABCout_RNAseq_Serum/model1_ENSMUSG00000089774_Serum_RNAseq_0.2_SMC_ratiotoproduction.txt"
counter = [0]
abs=[0]
final_abc =Array{Float64,2}(1000,4)
a=[]
b=[]
c=[]
d=[]
open(pathtofile_lisi) do file
    for ln in eachline(file)
        if startswith(ln, "-----")
            counter[1]=1
            abs[1]+=1
            # println("$(length(ln)), $(ln)")
        elseif counter[1]>0
            counter[1]+=1
            if counter[1] >= 5
                # println(ln)
                if (abs[1]==5)
                    ln_arr =split(ln)
                    push!(a,parse(Float64,ln_arr[1]))
                    push!(b,parse(Float64,ln_arr[2]))
                    push!(c,parse(Float64,ln_arr[3]))
                    push!(d,parse(Float64,ln_arr[4]))
                end
            end
        end
    end
end

a=a[1:end-2]
b=b[1:end-2]
c=c[1:end-2]
d=d[1:end-2]
a1=a
b1=b
c1=c
d1=d



# pathtofile ="/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/ABCout_RNAseq_Serum/model1_ENSMUSG00000089774_Serum_RNAseq_0.2_SMC_ratiotoproduction.txt"
counter = [0]
abs=[0]
final_abc =Array{Float64,2}(1000,4)
a=[]
b=[]
c=[]
d=[]
open(pathtofile_katie) do file
    for ln in eachline(file)
        if startswith(ln, "-----")
            counter[1]=1
            abs[1]+=1
            # println("$(length(ln)), $(ln)")
        elseif counter[1]>0
            counter[1]+=1
            if counter[1] >= 5
                # println(ln)
                if (abs[1]==5)
                    ln_arr =split(ln)
                    push!(a,parse(Float64,ln_arr[1]))
                    push!(b,parse(Float64,ln_arr[2]))
                    push!(c,parse(Float64,ln_arr[3]))
                    push!(d,parse(Float64,ln_arr[4]))
                end
            end
        end
    end
end

a=a[1:end-2]
b=b[1:end-2]
c=c[1:end-2]
d=d[1:end-2]



#
# fig = PyPlot.figure()
# sns.distplot(a,bins=20, kde= true, hist=true,rug= false ,  color = "orchid", label = "lisi",axlabel = "c act")
# sns.distplot(a1,bins=20, kde= true, hist=true,rug= false ,  color = "red", label = "katie",axlabel = "c act")
# PyPlot.savefig(string(folderpath,"c_act.pdf"))
#
# PyPlot.figure()
# sns.distplot(b,bins=20, kde= true, hist=true,rug= false ,  color = "orchid", label = "lisi",axlabel = "c de-act")
# sns.distplot(b1,bins=20, kde= true, hist=true,rug= false ,  color = "red", label = "katie",axlabel = "c de-act")
#
# PyPlot.savefig(string(folderpath,"c_de-act.pdf"))
#
# PyPlot.figure()
# sns.distplot(c,bins=20, kde= true, hist=true,rug= false ,  color = "orchid", label = "lisi",axlabel = "c deg")
# sns.distplot(c1,bins=20, kde= true, hist=true,rug= false ,  color = "red", label = "katie",axlabel = "c deg")
#
# PyPlot.savefig(string(folderpath,"c_deg.pdf"))
#
# PyPlot.figure()
# sns.distplot(d,bins=20, kde= true, hist=true,rug= false ,  color = "orchid", label = "lisi",axlabel = "scaling factor")
# sns.distplot(d1,bins=20, kde= true, hist=true,rug= false ,  color = "red", label = "katie",axlabel = "c act")
#
# PyPlot.savefig(string(folderpath,"scaling_factor.pdf"))



fig, ax = PyPlot.subplots(figsize=(20,15), ncols=3, nrows=4)
left   =  0.125  # the left side of the subplots of the figure
right  =  0.9    # the right side of the subplots of the figure
bottom =  0.1    # the bottom of the subplots of the figure
top    =  0.9    # the top of the subplots of the figure
wspace =  .5     # the amount of width reserved for blank space between subplots
hspace =  1.1    # the amount of height reserved for white space between subplots

# This function actually adjusts the sub plots using the above paramters
PyPlot.subplots_adjust(
    left    =  left,
    bottom  =  bottom,
    right   =  right,
    top     =  top,
    wspace  =  wspace,
    hspace  =  hspace
)
# The amount of space above titles
y_title_margin = 1.2
PyPlot.suptitle("Original vs Normalized vs Standardized", y = 1.09, fontsize=20)

### Bathrooms

sns.distplot(a,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c act",ax=ax[1,1])
sns.distplot(a,bins=20, kde= true, hist=false,rug= false ,  color = "orchid", label = "lisi",axlabel = "c act",ax=ax[1,2])
sns.distplot(a1,bins=20, kde= true, hist=false,rug= false ,  color = "red", label = "katie",axlabel = "c act",ax=ax[1,2])
sns.distplot(a1,bins=20, kde= true, hist=true,rug= true ,  color = "red", label = "katie",axlabel = "c act",ax=ax[1,3])

sns.distplot(b,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c de-act",ax=ax[2,1])
sns.distplot(b,bins=20, kde= true, hist=false,rug= false ,  color = "orchid", label = "lisi",axlabel = "c de-act",ax=ax[2,2])
sns.distplot(b1,bins=20, kde= true, hist=false,rug= false ,  color = "red", label = "katie",axlabel = "c de-act",ax=ax[2,2])
sns.distplot(b1,bins=20, kde= true, hist=true,rug= true ,  color = "red", label = "katie",axlabel = "c de-act",ax=ax[2,3])

sns.distplot(c,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c deg",ax=ax[3,1])
sns.distplot(c,bins=20, kde= true, hist=false,rug= false ,  color = "orchid", label = "lisi",axlabel = "c deg",ax=ax[3,2])
sns.distplot(c1,bins=20, kde= true, hist=false,rug= false ,  color = "red", label = "katie",axlabel = "c deg",ax=ax[3,2])
sns.distplot(c1,bins=20, kde= true, hist=true,rug= true ,  color = "red", label = "katie",axlabel = "c deg",ax=ax[3,3])

sns.distplot(d,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "scaling factor",ax=ax[4,1])
sns.distplot(d,bins=20, kde= true, hist=false,rug= false ,  color = "orchid", label = "lisi",axlabel = "scaling factor",ax=ax[4,2])
sns.distplot(d1,bins=20, kde= true, hist=false,rug= false ,  color = "red", label = "katie",axlabel = "scaling factor",ax=ax[4,2])
sns.distplot(d1,bins=20, kde= true, hist=true,rug= true ,  color = "red", label = "katie",axlabel = "scaling factor",ax=ax[4,3])
PyPlot.savefig(pathtoplot)
