using gene_bursting, Plots
sol = generate_single_simulation_m1([1.,1.,1.,1.],(0.,10.),base_model)
out = generate_single_simulation_samples_m1(
            [1.,1.,1.,1.],
            10,
            data_scaling_RNAseq)
plot(out)
