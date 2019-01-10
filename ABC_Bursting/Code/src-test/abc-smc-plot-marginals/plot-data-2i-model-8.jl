using PyCall, PyPlot; @pyimport seaborn as sns
#example call : julia plot-data-2i.jl ENSMUSG00000029472
# name = "ENSMUSG00000006398"
# pathtofile_lisi = string(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_8_2i_easy/model_8_easy_",name),"_2i_0.8_100.txt")
# pathtoplot = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/model-8-easy/",name)
# namers=["ENSMUSG00000006398", "ENSMUSG00000027715", "ENSMUSG00000029177", "ENSMUSG00000029472", "ENSMUSG00000030867", "ENSMUSG00000041431"]
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
    pathtofile_lisi = string(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/lisi_high_var_2i/model1_",name),"_2i_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt")
    pathtofile_katie = string(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/ABCout_RNAseq_transformed_keygenes_2_2i/model1_",name),"_2i_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt")
    pathtoplot= ARGS[2]
elseif length(ARGS) == 1
    print("you have passed 1 argument. i am thinking it is gene name:")
    name = ARGS[1]
    pathtofile_lisi = string(string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/model_8_2i_easy/model_8_easy_",name),"_2i_0.8_100.txt")
    pathtoplot = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/model-8-easy/",name)
    # println(name)
    # println(pathtofile_lisi)
    # println(pathtoplot)

end



counter_smth = [0]
counter_abcsmcs=[0]
count_particles = [0]
#enter number of accepted particles and number of unknown parameters
final_abc =Array{Float64,2}(100,5)
a=[]
b=[]
c=[]
d=[]
e=[]
open(pathtofile_lisi) do file
    for ln in eachline(file)
        if startswith(ln, "-----")
            counter_smth[1]=1
            counter_abcsmcs[1]+=1
            # println("$(length(ln)), $(ln)")
        elseif counter_smth[1]>0
            counter_smth[1]+=1
            if counter_smth[1] >= 5
                if (counter_abcsmcs[1]==3)
                    count_particles[1]+=1
                    if count_particles[1]<=100

                        ln_arr =split(ln)
                        push!(a,parse(Float64,ln_arr[1]))
                        push!(b,parse(Float64,ln_arr[2]))
                        push!(c,parse(Float64,ln_arr[3]))
                        push!(d,parse(Float64,ln_arr[4]))
                        push!(e,parse(Float64,ln_arr[5]))
                    end
                end
            end
        end
    end
end
a=a[1:end-2]
b=b[1:end-2]
c=c[1:end-2]
d=d[1:end-2]
e=e[1:end-2]
PyPlot.ion()
fig, ax = PyPlot.subplots(figsize=(6,15),ncols=1, nrows=5)
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

sns.distplot(a,bins=20,  kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c act",ax=ax[1,1])
sns.distplot([-4,4],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",axlabel = "c act",ax=ax[1,1])
PyPlot.suptitle("Original vs Normalized vs Standardized", y = 1.09, fontsize=20)
sns.distplot(b,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c de-act",ax=ax[2,1])
sns.distplot([-4,4],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",axlabel = "c de-act",ax=ax[2,1])

sns.distplot(c,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c deg",ax=ax[3,1])
sns.distplot([-4,4],bins=20,  kde= false, hist=true,rug= false ,  color = "white",  label = "lisi",axlabel = "c deg",ax=ax[3,1])

sns.distplot(d,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "feedback",ax=ax[4,1])
sns.distplot([-4,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white",  label = "lisi",axlabel = "feedback",ax=ax[4,1])

sns.distplot(e,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "scaling factor",ax=ax[5,1])
sns.distplot([-1,6],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",axlabel = "scaling factor",ax=ax[5,1])


PyPlot.savefig(pathtoplot)


Distributions.Uniform(-3.0,3.0),
    Distributions.Uniform(-3.0,3.0),
    Distributions.Uniform(-3.0,3.0),
    Distributions.Uniform(-3.0,0.5),
    Distributions.Uniform(0.0,5.0)
