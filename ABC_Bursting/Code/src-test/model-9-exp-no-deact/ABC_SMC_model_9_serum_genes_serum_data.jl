@everywhere using ABC_Bursting
# CHANGE: outname and epsilons
println("Start include.")
@everywhere include("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Code/external/ABaCus.jl")
println("Starting")
gene_names = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/gene_names/genes_high_var_serum.txt",'\t','\n')
# gene_names = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/gene_names/genes_high_var_serum.txt",'\t','\n')

gene_names = gene_names[1:end,:,:]
data = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_serum_filtereddata.txt")
# data = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_serum_filtereddata.txt")
names = data[:,1]
gene_names = intersect(data[2:end,1],gene_names)
println("Ready")
selectionstart = parse(Int64,ARGS[1])
selectionend = parse(Int64,ARGS[2])
outfile =  ARGS[3]
@sync @parallel for i in gene_names[selectionstart:selectionend]
    test_data = Get_single_gene_RNAseq_data(i,data)
    outfile_name = "$(outfile)/model_9_serum_genes_serum_data_$(i)_small_priors_02.txt"
    n_simulations = length(test_data)
    priors = [Distributions.Uniform(-6.0,0.0),
        Distributions.Uniform(-6.0,0.0),
        Distributions.Uniform(-5.0,2.0),#b
        Distributions.Uniform(0.0,5.0)] #k
    n_params = 4 #number of parameters to be inferred
    n_particles = 100 #number of particles accepted
    threshold = [0.8,0.4,0.2]

    println("Parameters set for $i")
    SMC_input = ABaCus.ABCSMCInput(
                                        n_params,
                                        n_particles,
                                        threshold,
                                        priors,
                                        kolmogorov_smirnov_distance,
                                        generate_single_simulation_samples_m9_no_deact
                                        )

    println("starting ABC for $i")
    abacus_output_para_est_one_gene = @time ABaCus.ABCSMC(SMC_input,
                                test_data,
                                n_simulations,
                                data_scaling_RNAseq
                                )


    # PyPlot.figure()
    # M = [abacus_output_para_est_one_gene.population[1]';abacus_output_para_est_one_gene.population[2]';abacus_output_para_est_one_gene.population[3]']
    # #M = [result.population[1]';result.population[2]']
    # #figure1 = corrplot(M, label = ["x$i" for i=1:4])
    # #Plots.savefig(figure1,"myfigure1.png")
    # for j in 1:4
    #    pp = M[:,j]'
    #    PyPlot.subplot(140+j)
    #    PyPlot.plt[:hist](pp,5)
    #
    # end
    println("Writing $i to file")
    ABaCus.write(outfile_name, abacus_output_para_est_one_gene)
end
