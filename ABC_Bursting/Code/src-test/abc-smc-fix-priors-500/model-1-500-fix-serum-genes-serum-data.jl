################################################# running for serum ######################################################
@everywhere using ABC_Bursting
println("Start include.")
@everywhere include("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Code/external/ABaCus.jl")
println("Starting")
gene_names = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/gene_names/genes_high_var_serum.txt",'\t','\n')
gene_names = gene_names[1:end,:,:]
data = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_serum_filtereddata.txt")
names = data[:,1]
gene_names = intersect(data[2:end,1],gene_names)
println("Ready")
selectionstart = parse(Int64,ARGS[1])
selectionend = parse(Int64,ARGS[2])
outfile =  ARGS[3]
@sync @parallel for i in gene_names[selectionstart:selectionend]
    test_data = Get_single_gene_RNAseq_data(i,data)
    outfile_name = "$(outfile)/model1_$(i)_0.1_500_250000_serum_genes_serum_data.txt"
    n_simulations = length(test_data)
    priors = [Distributions.Uniform(-6.0,0.0),
        Distributions.Uniform(-6.0,0.0),
        Distributions.Uniform(-6.0,0.0)]
    n_params = 3
    n_particles = 500
    threshold = [0.4,0.2,0.1]
    println("Parameters set for $i")
    SMC_input = ABaCus.ABCSMCInput(
                                        n_params,
                                        n_particles,
                                        threshold,
                                        priors,
                                        kolmogorov_smirnov_distance,
                                        generate_single_simulation_samples_m1_250000
                                        )

    println("starting ABC for $i")
    abacus_output_para_est_one_gene = @time ABaCus.ABCSMC(SMC_input,
                                test_data,
                                n_simulations,
                                data_scaling_RNAseq
                                )

    println("Writing $i to file")
    ABaCus.write(outfile_name, abacus_output_para_est_one_gene)
end
