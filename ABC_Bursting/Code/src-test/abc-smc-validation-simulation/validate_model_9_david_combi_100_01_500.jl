@everywhere using ABC_Bursting
# CHANGE: outname and epsilons
println("Start include.")
@everywhere include("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Code/external/ABaCus.jl")
println("Starting")
gene_names = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/gene_names/genes_high_var_2i.txt",'\t','\n')
# gene_names = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/gene_names/genes_high_var_serum.txt",'\t','\n')

gene_names = gene_names[1:end,:,:]
data = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_2i_filtereddata.txt")
# data = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_serum_filtereddata.txt")
names = data[:,1]
gene_names = intersect(data[2:end,1],gene_names)

# Completed_SCPopnSerum_serum = Get_completed_genes_list_RNAseq(gene_names,"/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/ABCout_RNAseq_transformed_keygenes_2_serum/model1_","_serum_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt")
# Completed_SCPopnSerum_serum = Completed_SCPopnSerum_serum[1]
# SCPopnSerum_serum = intersect(gene_names,Completed_SCPopnSerum_serum)
#Only runs gene names which are present in the data file
# gene_names = setdiff(gene_names,SCPopnSerum_serum)

println("Ready")

#Takes input for selection of which genes are to be run
#This allows for separation of the list to run across different CPUs)
#----------------------------------------------------------------
selectionstart = parse(Int64,"1")
selectionend = parse(Int64,"1")
outfile = "/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/ABC_output/RNAseq/validation_simulation/model_9"


# selectionstart = parse(Int64,"1")
# selectionend = parse(Int64,"2")
# outfile = "tester"

#Run ABC for selected genes and the experiemental data
#----------------------------------------------------------------
@sync @parallel for i in gene_names[selectionstart:selectionend]
    test_data = generate_single_simulation_samples_m9([-3.,-2.3,-4.,-2.3,2],500, data_scaling_RNAseq)

    #Define output file
    #----------------------------------------------------------------
    outfile_name = "$(outfile)/validate_model_9_david_100_01_500.txt"
    n_simulations = length(test_data)
    #Define parameters for the ABC run
    #----------------------------------------------------------------
    #Define priors
    priors = [Distributions.Uniform(-6.0,0.0),
        Distributions.Uniform(-6.0,0.0),
        Distributions.Uniform(-6.0,0.0),
        Distributions.Uniform(-20.0,10.0),
        Distributions.Uniform(-20.0,20.0)] #changes from 0.0,5.0
    #Define parameters for SMC
    n_params = 5 #number of parameters to be inferred
    n_particles = 100 #number of particles accepted
    threshold = [0.8,0.4,0.1]

    println("Parameters set for $i")
    SMC_input = ABaCus.ABCSMCInput(
                                        n_params,
                                        n_particles,
                                        threshold,
                                        priors,
                                        kolmogorov_smirnov_distance,
                                        generate_single_simulation_samples_m9
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
