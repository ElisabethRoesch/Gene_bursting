using ABC_Bursting, DifferentialEquations
data = readdlm("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/Data/data_files/RNAseq_2i_filtereddata.txt")
res_i = Get_single_gene_RNAseq_data("ENSMUSG00000029472",data)

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

res_j = generate_single_simulation_samples_m1(
            [-2.8,-2.7,-2.5,4],
            295,
            data_scaling_RNAseq
            )


using Plots
histogram(res_i)
histogram(res_j)

test_data = Get_single_gene_RNAseq_data(i,data)

outfile_name = "$(outfile)/model1_cellcycle_$(i)_2i_0.1_100.txt"
n_simulations = length(test_data)

priors = [Distributions.Uniform(-3.0,3.0),
    Distributions.Uniform(-3.0,3.0),
    Distributions.Uniform(-3.0,0.5),
    Distributions.Uniform(0.0,5.0)]

n_params = 4
n_particles = 100
threshold = [0.8,0.4,0.1]

println("Parameters set for $i")
SMC_input = ABaCus.ABCSMCInput(
                                    n_params,
                                    n_particles,
                                    threshold,
                                    priors,
                                    kolmogorov_smirnov_distance,
                                    cell_cycle_simulation_generate_single_simulation_samples_m1
                                    )

println("starting ABC for $i")
abacus_output_para_est_one_gene = @time ABaCus.ABCSMC(SMC_input,
                            test_data,
                            n_simulations,
                            data_scaling_RNAseq
                            )

println("Writing $i to file")
ABaCus.write(outfile_name, abacus_output_para_est_one_gene)
