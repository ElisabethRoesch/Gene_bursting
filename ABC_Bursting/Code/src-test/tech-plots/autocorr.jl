#testing for basic understanding

using ABC_Bursting, DifferentialEquations

base_model = @reaction_network CC begin
          c1, G0 --> G1
          c2, G1 --> G0
          c3, G1 --> G1 + P
          c4, P --> 0
end c1 c2 c3 c4
function data_scaling_RNAseq(outcome::Array{Int64}, parameter::Float64)
    output = log10.((outcome[3] * parameter).+1)
    return output
end
#act,de-act,deg
max_marg = [1.,-3,-3,-1.8,5.]
parameter_input = max_marg
n_samples=1000




c1 = 10^parameter_input[1]  #Activation
c2 = 10^parameter_input[2]  #Deactivation
c3 = 1.0                    #Expression
c4 = 10^parameter_input[3]  #Degradation
k = 10^parameter_input[4]     #Feedback
c6 = parameter_input[5]     #Scale factor
#Define starting state at stationary distribution
G = 2.0
P = (sqrt((4*c1*c2*c3*c4*G*k)+(c1*c4*k+c2*c4*k)^2)-c1*c4*k-c2*c4*k)/(2*c2*c4)
G1 = round(Int64,(
          (c1*G*k)/((c1*k)+(c2*(k+P)))
          ))
P = round(Int64,(P))
if P < 0
          P = 0
else
          P = P
end
if G1 < 0
          G1 = 0
else
          G1 = G1
end
G0 = 2 - G1
starting_state = [G0,G1,P]
time = 40000.0
time_range= (0.0,time)
#Define model
#Set rate parameters

rate1 = function (u,p,t)
           return (c1)*u[1]*(1/(1+(u[3]/k)))
end

affect! = function (integrator)
          integrator.u[1] -= 1
          integrator.u[2] += 1
end
jump1 = DifferentialEquations.ConstantRateJump(rate1,affect!)

rate2 = function (u,p,t)
          return (c2)*u[2]
end

affect! = function (integrator)
          integrator.u[1] += 1
          integrator.u[2] -= 1
end
jump2 = DifferentialEquations.ConstantRateJump(rate2,affect!)

rate3 = function (u,p,t)
          return (c3)*u[2]
end
affect! = function (integrator)
          integrator.u[2] += 0
          integrator.u[3] += 1
end
jump3 = DifferentialEquations.ConstantRateJump(rate3,affect!)

rate4 = function (u,p,t)
          return (c4)*u[3]
end

affect! = function (integrator)
          integrator.u[3] -= 1
end
jump4 = DifferentialEquations.ConstantRateJump(rate4,affect!)

prob = DifferentialEquations.DiscreteProblem(starting_state,time_range)
jump_prob = DifferentialEquations.JumpProblem(prob,DifferentialEquations.Direct(),jump1,jump2,jump3,jump4)

#Simulate multiple cells and take end point
distribution_data = Array{Float64}(n_samples)
sol = DifferentialEquations.solve(jump_prob,FunctionMap())
max_time = sol.t[end]
for i in 1:n_samples
          outcome = round.(Int64,(sol(i*(max_time/(n_samples+1)))))
          distribution_data[i,1] = data_scaling_RNAseq(outcome, c6)
end

r=range(1,n_samples)
using Plots
plot(sol,xlab="time",ylab="product abundance", title="slow: sol of jumpproblem")
savefig("idea10.pdf")
plot(r,distribution_data,title ="slow: product abundance of 100 cells/timepoints converted by scaling factor",xlab="cell/timepoint",ylab="product abundance",label =["P"])
title =string("simulated product abundance of ", n_samples," cells")
histogram(distribution_data,label=["model 8 max. marginal product abundance"],title=title, xlab="product abundance", ylab = "frequency",nbins=20)
savefig("ideahist10.pdf")


parameter_input = [-2.3,-1.1,-1.7]
parameter_input = [-2.3,-2.9,-2.2] #auto egal

w=log10.([0.001, 0.0001, 0.04, 1000])

function data_scaling_RNAseq(outcome::Array{Int64}, parameter::Float64)
    output = log10.((outcome[3] * parameter).+1)
    return output
end
parameter_input = [-2.3,-1.1,-1.7]
as=[]
for i in 1:10
          push!(as,generate_single_simulation_samples_m1_50000( parameter_input,125,  data_scaling_RNAseq))
end
b=generate_single_simulation_samples_m1_100000( parameter_input,250,  data_scaling_RNAseq)
c=generate_single_simulation_samples_m1_125000( parameter_input,250,  data_scaling_RNAseq)
d=generate_single_simulation_samples_m1_200000( parameter_input,250,  data_scaling_RNAseq)
e=generate_single_simulation_samples_m1_250000( parameter_input,125,  data_scaling_RNAseq)
f=generate_single_simulation_samples_m1_500000( parameter_input,125,  data_scaling_RNAseq)



using PyCall, PyPlot;
@pyimport seaborn as sns
fig, ax = PyPlot.subplots()
sns.distplot(a,bins=15, kde= true, hist=false,rug= false , label = "125 of 50000")
# sns.distplot(b,bins=15, kde= true, hist=false,rug= false ,  label = "100000")
# sns.distplot(c,bins=15, kde= true, hist=false,rug= false ,   label = "125000")
# sns.distplot(d,bins=15, kde= true, hist=false,rug= false ,   label = "200000")
# sns.distplot(e,bins=15, kde= true, hist=false,rug= false ,   label = "125 of 250000")
# sns.distplot(f,bins=15, kde= true, hist=false,rug= false ,   label = "125 of 500000")
ax[:set_ylabel]("Density",fontsize=16)
ax[:set_xlabel]("Product abundance",fontsize=16)
ax[:tick_params](labelsize=12)
ax[:legend](fontsize=12)
xs=[]
ys=[]
zs=[]
mini = min(min(as...)...)
max = max(max(as...)...)
for i in 1:10
          x,y,z = PyPlot.plt[:hist](a,20,density=true,)
          push!(xs,x)
          push!(ys,y)
          push!(zs,z)
end
xs


using PyPlot
using ABC_Bursting
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport pandas as pd
using PyCall, PyPlot; @pyimport numpy as np
# ENSMUSG00000025056 ENSMUSG00000048402 ENSMUSG00000021835 ENSMUSG00000048402
# genes = [ "ENSMUSG00000021835","ENSMUSG00000025056",
genes =["ENSMUSG00000048402"]
# genes = ["ENSMUSG00000051176","ENSMUSG00000020717","ENSMUSG00000048402"]
sampled_parameter_combinations_serum = []
n_samples = 500
n_tester = 100
binsize = 20
#read experimental data_2i
all_genes_min_max_serum = []
all_genes_experimental_results_serum = []
all_genes_all_simulation_results_serum = []
all_genes_samples_serum = []
data_serum = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_2i_filtereddata.txt")
for name in genes
      path_to_file_serum =  string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/explore_distr/model_1/2i/model1_",name, "_2i_0.1_fix_500_250000.txt")
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
for i in 1:length(genes)
    gene = genes[i]
    #read and store experimental data_2i for gene "gene"
    one_gene_experimental_result_serum = Get_single_gene_RNAseq_data(gene,data_serum)
    push!(all_genes_experimental_results_serum,one_gene_experimental_result_serum)
    #simulate and do histogram of samples
    one_gene_all_simulation_results_serum= []
    one_gene_bins_serum = Dict()
    one_gene_bins_serum["bins"] =  range(1,binsize)
    for sim in 1:length(sampled_parameter_combinations_serum[i])
        one_gene_one_simulation_result_serum= generate_single_simulation_samples_m1_50000(  [-2.3,-1.1,-1.7],125,  data_scaling_RNAseq)
        push!(one_gene_all_simulation_results_serum, one_gene_one_simulation_result_serum)
    end
    min_ele_serum = min(min(one_gene_all_simulation_results_serum...)...)
    max_ele_serum = max(max(one_gene_all_simulation_results_serum...)...)
    push!(all_genes_min_max_serum, [min_ele_serum, max_ele_serum])
    for  sim in 1:length(sampled_parameter_combinations_serum[i])
        bins_serum, indata_serum =np.histogram(one_gene_all_simulation_results_serum[sim], density = true, bins = binsize, range = (min_ele_serum, max_ele_serum))
        one_gene_bins_serum["sample_$(sim)"] = bins_serum
    end
    push!(all_genes_samples_serum,one_gene_bins_serum)
    push!(all_genes_all_simulation_results_serum,one_gene_all_simulation_results_serum)
end

for gene in 1:length(genes)
    data_serum = all_genes_samples_serum[gene]

    in_para_serum =[]

    colnames_serum = []

    for key in keys(all_genes_samples_serum[gene])
        colname=(string(key))
        push!(colnames_serum,colname)
    end

    df_serum = pd.DataFrame(data_serum, columns = colnames_serum)


    for i in 1:length(colnames_serum)
        push!(in_para_serum,df_serum[colnames_serum[i]])
    end

    fig,ax = PyPlot.subplots()

      sns.tsplot(in_para_serum,color="salmon", time=np.linspace(all_genes_min_max_serum[gene][1], all_genes_min_max_serum[gene][2], binsize),  ci=75,ax=ax)
      # sns.distplot(all_genes_experimental_results_serum[gene],norm_hist=true, bins=binsize, kde= true, hist=false,rug= false ,  color = "red", label = "experimental 2i")
      ax[:set_ylabel]("Frequency",fontsize=16)
      ax[:set_xlabel]("Product abundance",fontsize=16)
      name = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/explore-distr/2i_rightpeak/", genes[gene],"_n_tester",string(n_tester),"_binsize",string(binsize),"2.pdf")
     # savefig(name)
end
