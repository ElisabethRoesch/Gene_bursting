using PyPlot
using ABC_Bursting
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport pandas as pd
using PyCall, PyPlot; @pyimport numpy as np

genes = [ "ENSMUSG00000012396"]
data_2i = []
data_serum = []
n_samples = 100
tau = 1000
tmax = tau*n_samples

fileending_2i = string("_2i_0.1_fix_",n_samples,"_",tmax)
fileending_serum = string("_serum_0.1_fix_",n_samples,"_",tmax)
for name in genes
    path_to_plot = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/model-1-nanog/500/serum/",name,fileending_2i,fileending_serum)
    path_to_file_2i = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/",tau,"/2i/model1_",name,fileending_2i,".txt")
    path_to_file_serum = string("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/testing_priors/model_1/nanog/",tau,"/serum/model1_",name,fileending_serum,".txt")
    data_file_2i=read_one_m1(path_to_file_2i,  3, n_samples,3)
    data_file_serum=read_one_m1(path_to_file_serum,  3, n_samples,3)
    push!(data_2i, data_file_2i)
    push!(data_serum, data_file_serum)
    cum_posterior_2i = Vector(n_samples)
    cum_posterior_serum = Vector(n_samples)
    cumsum!(cum_posterior_2i, data_file_2i[6])
    cumsum!(cum_posterior_serum, data_file_serum[6])
    sum_2i = string(sum(data_2i[1][5]))
    sum_serum = string(sum(data_serum[1][5]))

    fig,ax = PyPlot.subplots()
        PyPlot.plot( cum_posterior_2i,".",color="salmon")
        ax[:set_ylabel]("cdf")
        ax[:set_xlabel]("Accepted parameters")
    savefig(string(path_to_plot,"cdf_2i.png"))

    fig,ax = PyPlot.subplots()
        PyPlot.plot( cum_posterior_serum,".",color="cornflowerblue")
        ax[:set_ylabel]("cdf")
        ax[:set_xlabel]("Accepted parameters")
    savefig(string(path_to_plot,"cdf_serum.png"))

    fig,ax = PyPlot.subplots()
        sns.distplot(data_serum[1][5],color="cornflowerblue")
        ax[:set_ylabel]("Frequency")
        ax[:set_xlabel]("Accepted distance")
    savefig(string(path_to_plot,"distances_2i.png"))

    fig,ax = PyPlot.subplots()
        sns.distplot(data_2i[1][5],color="salmon")
        ax[:set_ylabel]("Frequency")
        ax[:set_xlabel]("Accepted distance")
    savefig(string(path_to_plot,"distances_serum.png"))

end
