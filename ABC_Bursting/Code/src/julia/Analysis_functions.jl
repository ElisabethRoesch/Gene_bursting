 import StatsBase
import Gadfly
include("./Code/Gillespie_model_1_fixedprod.jl")
include("./Code/ABC_functions.jl")
include("./Code/ABaCus.jl")

#Samples from the joint posteriors
function Posterior_sampling(gene_name,file_prefix,file_suffix,n_samples)
  i = gene_name
  test_output = ABaCus.read_smc_output("$(file_prefix)$(i)$(file_suffix)")
  params = test_output.population[5]
  weights = test_output.weights[5] #5 - set by the number of thresholds used in the SMC

  store = Array{Float64}(n_samples,4) #4 - set by the number of parameters being inferred
  for i in 1:n_samples
    index = StatsBase.sample(collect(1:1000), weights) #1000 - set by the number of particles accepted in ABC
    store[i,:] = params[:,index]
  end
  return store
end

#Plots the product levels for an individual run
function plot_single_simulation(sol,colour_pick)
    cols = length(sol.u)
    rows = length(sol.u[1])
    sol_data = Array{Float64}(cols,rows)

    for i in 1:cols
        for j in 1:length(sol.u[1])
            sol_data[i,j] = sol.u[i][j]
        end
    end
    Gadfly.plot(x= sol.t ,y = sol_data[:,3], Gadfly.Geom.line, Gadfly.Theme(default_color=colour_pick),
    Gadfly.Guide.xlabel("Time (AU)"), Gadfly.Guide.ylabel("Product Abundance (molecules)"))
end

#Plots a density of the simulated data
function plot_simulated_distribution_dens(data)
    Gadfly.plot(x=data , Gadfly.Geom.density())
end

#Plots comparison between two simulations e.g. compare data to a simulation
function plot_simulated_distribution_comparison(data1, data2)
    Gadfly.plot(Gadfly.layer(x=data1 , Gadfly.Geom.density(),Gadfly.Theme(default_color="red")),
    Gadfly.layer(x=data2 , Gadfly.Geom.density(),Gadfly.Theme(default_color="blue")),
    Gadfly.Guide.manual_color_key("Data", ["Data1", "Data2"], ["red", "blue"]),
           Gadfly.Guide.xlabel("RNAser output"),
           Gadfly.Guide.ylabel("fraction of cells"))
end

#Plot parameter posteriors
function plot_parameters_4(p1::Array{Float64},p2::Array{Float64},p3::Array{Float64},p4::Array{Float64})
    P = Array{Float64}(length(p1),4)
    P[:,1] = p1
    P[:,2] = p2
    P[:,3] = p3
    P[:,4] = p4
    StatPlots.corrplot(P)
end

function plot_output_parameters(gene,file_location,particle_num)

            test_output = ABaCus.read_smc_output(file_location)

            p1_posterior = Array{Float64}(particle_num)
            p2_posterior = Array{Float64}(particle_num)
            p3_posterior = Array{Float64}(particle_num)
            p4_posterior = Array{Float64}(particle_num)

            x = 5
            for j in 1:particle_num
                p1_posterior[j] = test_output.population[x][:,j][1]
                p2_posterior[j] = test_output.population[x][:,j][2]
                p3_posterior[j] = test_output.population[x][:,j][3]
                p4_posterior[j] = test_output.population[x][:,j][4]
            end
            return plot_parameters_4(p1_posterior,p2_posterior,p3_posterior,p4_posterior)
end


function Plot_multi_plot_serun_2i(posterior_sample_serum,posterior_sample_2i)

	data_serum = readdlm("./Data/data_files/RNAseq_serum_filtereddata.txt")
	data_2i = readdlm("./Data/data_files/RNAseq_2i_filtereddata.txt")

	nanog_data_serum = Get_single_gene_RNAseq_data(gene_name,data_serum)
	nanog_data_2i = Get_single_gene_RNAseq_data(gene_name,data_2i)

	Act_test_serum = posterior_sample_serum[:,1]-posterior_sample_serum[:,3]
	Deact_test_serum = posterior_sample_serum[:,2]-posterior_sample_serum[:,3]

	Act_test_2i = posterior_sample_2i[:,1]-posterior_sample_2i[:,3]
	Deact_test_2i= posterior_sample_2i[:,2]-posterior_sample_2i[:,3]


	nanog_params_act = Gadfly.plot(Gadfly.layer(x=posterior_sample_serum[:,1], Gadfly.Geom.density,Gadfly.Theme(default_color="blue")),
				Gadfly.layer(x=posterior_sample_2i[:,1], Gadfly.Geom.density,Gadfly.Theme(default_color="red")),
		Gadfly.Guide.xlabel("Activation rate (log)"),
	        Gadfly.Guide.ylabel("Density"),
		Gadfly.Theme(key_position = :none),
		Gadfly.Coord.Cartesian(xmin=-3,xmax=3))

	nanog_params_deact = Gadfly.plot(Gadfly.layer(x=posterior_sample_serum[:,2], Gadfly.Geom.density,Gadfly.Theme(default_color="blue")),
				Gadfly.layer(x=posterior_sample_2i[:,2], Gadfly.Geom.density,Gadfly.Theme(default_color="red")),
		Gadfly.Guide.xlabel("Deactivation rate (log)"),
	        Gadfly.Guide.ylabel("Density"),
		Gadfly.Theme(key_position = :none),
		Gadfly.Coord.Cartesian(xmin=-3,xmax=3))

	nanog_params_deg = Gadfly.plot(Gadfly.layer(x=posterior_sample_serum[:,3], Gadfly.Geom.density,Gadfly.Theme(default_color="blue")),
				Gadfly.layer(x=posterior_sample_2i[:,3], Gadfly.Geom.density,Gadfly.Theme(default_color="red")),
		Gadfly.Guide.xlabel("Degradation rate (log)"),
	        Gadfly.Guide.ylabel("Density"),
		Gadfly.Theme(key_position = :none),
		Gadfly.Coord.Cartesian(xmin=-3,xmax=0.5))

	nanog_params_scale = Gadfly.plot(Gadfly.layer(x=posterior_sample_serum[:,4], Gadfly.Geom.density,Gadfly.Theme(default_color="blue")),
				Gadfly.layer(x=posterior_sample_2i[:,4], Gadfly.Geom.density,Gadfly.Theme(default_color="red")),
		Gadfly.Guide.xlabel("Scale factor"),
	        Gadfly.Guide.ylabel("Density"),
		Gadfly.Theme(key_position = :none),
		Gadfly.Coord.Cartesian(xmin=0,xmax=5))


	sim_dist_store_serum = Array{Float64}(250,100)
	for i in 1:100
		dist_serum = generate_single_simulation_samples_m1(posterior_sample_serum[i,:],250,data_scaling_RNAseq,cell_cycle_model)
		sim_dist_store_serum[:,i] = dist_serum
	end

	sim_dist_store_2i = Array{Float64}(250,100)
	for i in 1:100
		dist_2i = generate_single_simulation_samples_m1(posterior_sample_2i[i,:],250,data_scaling_RNAseq,cell_cycle_model)
		sim_dist_store_2i[:,i] = dist_2i
	end

	nanog_simdata_plot = Gadfly.plot(Gadfly.layer(x=nanog_data_serum, Gadfly.Geom.density,Gadfly.Theme(default_color="blue")),
			Gadfly.layer(x=sim_dist_store_serum[:,1], Gadfly.Geom.density,Gadfly.Theme(default_color="green")),
			Gadfly.layer(x=sim_dist_store_serum[:,2], Gadfly.Geom.density,Gadfly.Theme(default_color="green")),
			Gadfly.layer(x=sim_dist_store_serum[:,3], Gadfly.Geom.density,Gadfly.Theme(default_color="green")),
			Gadfly.layer(x=sim_dist_store_serum[:,4], Gadfly.Geom.density,Gadfly.Theme(default_color="green")),
			Gadfly.layer(x=sim_dist_store_serum[:,5], Gadfly.Geom.density,Gadfly.Theme(default_color="green")),
			Gadfly.layer(x=nanog_data_2i, Gadfly.Geom.density,Gadfly.Theme(default_color="red")),
			Gadfly.layer(x=sim_dist_store_2i[:,1], Gadfly.Geom.density,Gadfly.Theme(default_color="orange")),
			Gadfly.layer(x=sim_dist_store_2i[:,2], Gadfly.Geom.density,Gadfly.Theme(default_color="orange")),
			Gadfly.layer(x=sim_dist_store_2i[:,3], Gadfly.Geom.density,Gadfly.Theme(default_color="orange")),
			Gadfly.layer(x=sim_dist_store_2i[:,4], Gadfly.Geom.density,Gadfly.Theme(default_color="orange")),
			Gadfly.layer(x=sim_dist_store_2i[:,5], Gadfly.Geom.density,Gadfly.Theme(default_color="orange")),
		Gadfly.Guide.xlabel("RNAseq value"),
	        Gadfly.Guide.ylabel("density"),
		Gadfly.Theme(key_position = :none))

	Horizontalx = [-3,3]
	Horizontaly = [0,0]
	Verticalx = [0,0]
	Verticaly = [-3,0]

	nanog_param_norm = Gadfly.plot(Gadfly.layer(x=Deact_test_serum, y=Act_test_serum, Gadfly.Geom.point,Gadfly.Theme(default_color="blue")),
				Gadfly.layer(x=Deact_test_2i, y=Act_test_2i ,Gadfly.Geom.point,Gadfly.Theme(default_color="red")),
				Gadfly.layer(x=Horizontalx, y=Horizontaly ,Gadfly.Geom.line,Gadfly.Theme(default_color="black")),
				Gadfly.layer(x=Verticalx, y=Verticaly ,Gadfly.Geom.line,Gadfly.Theme(default_color="black")),
		Gadfly.Guide.xlabel("Deact/Deg"),
	        Gadfly.Guide.ylabel("Act/Deg"),
		Gadfly.Theme(key_position = :none))



	left_stack = Gadfly.vstack(nanog_simdata_plot,nanog_param_norm)
	right_stack = Gadfly.vstack(nanog_params_act,nanog_params_deact,nanog_params_deg,nanog_params_scale)
	stack = Gadfly.hstack(left_stack,right_stack,)


	title = compose(context(0, 0, 1, 0.1),
                text(0.5, 1.0, names_dict[gene_name], hcenter, vbottom))
	p = vstack(title, compose(context(0, 0, 1, 0.9), stack))
    Gadfly.draw(Gadfly.SVG("./Plots/Selected_genes_ser2i_summary_plots/$(names_dict[gene_name])_posterior_summary.svg", 6Gadfly.inch, 6Gadfly.inch), p)
end


#Demonstration________________________________________________________________________

#plot single simulation through time
sol = generate_single_simulation_m1(posterior_sample_serum[1,:],(0.0,40000.0),cell_cycle_model)
plot_single_simulation(sol,"blue")

#plot sampled distribution
sim_samples = generate_single_simulation_samples_m1(posterior_sample_2i[1,:],250,data_scaling_RNAseq,cell_cycle_model)
plot_simulated_distribution_dens(sim_samples)

#plot sampled distribution vs real data
sim_samples = generate_single_simulation_samples_m1(posterior_sample_2i[1,:],250,data_scaling_RNAseq,cell_cycle_model)
data_2i = readdlm("./Data/data_files/RNAseq_2i_filtereddata.txt")
nanog_data_2i = Get_single_gene_RNAseq_data(gene_name,data_2i)
plot_simulated_distribution_comparison(nanog_data_2i, sim_samples)

#plot posteriors (big plotting function)
names = readdlm("./Data/data_files/Selected_genes_names.txt",'\t','\n') #This needs checking/fixing
names = names[:,1:2] #This needs checking/fixing
names_dict = Dict(names[i,1] => names[i,2] for i = 1:length(names[:,1]))

posterior_sample_serum = Posterior_sampling(gene_name,"./Data/ABC_output/RNAseq/ABCout_RNAseq_transformed_keygenes_2_serum/model1_","_Serum_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt",200)
posterior_sample_2i = Posterior_sampling(gene_name,"./Data/ABC_output/RNAseq/ABCout_RNAseq_transformed_keygenes_2_2i/model1_","_2i_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt",200)
Plot_multi_plot_serun_2i(posterior_sample_serum,posterior_sample_2i) #Saves plot to file defined in function



#_______________________________________________________________________________
#Plotting of the large plot for nanog found in the paper: nanog posteriors -  Plotting in R
#The data saved here is used to plot in R (Posterior_plotting.R)


#posterior plots
gene_name = "ENSMUSG00000012396"
posterior_sample_serum = Posterior_sampling(gene_name,"./Data/ABC_output/RNAseq/ABCout_RNAseq_transformed_keygenes_2_serum/model1_","_Serum_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt",100)
posterior_sample_2i = Posterior_sampling(gene_name,"./Data/ABC_output/RNAseq/ABCout_RNAseq_transformed_keygenes_2_2i/model1_","_2i_RNAseq_0.122_SMC_ratiotoproduction_transformed.txt",100)
writedlm("./Data/Analysis_output/nanog_serum_posterior_sample.txt", posterior_sample_serum, '\t')
writedlm("./Data/Analysis_output/nanog_2i_posterior_sample.txt", posterior_sample_2i, '\t')

#normalise parameters by degradation rate
Act_test_serum = posterior_sample_serum[:,1]-posterior_sample_serum[:,3]
Deact_test_serum = posterior_sample_serum[:,2]-posterior_sample_serum[:,3]
Act_test_2i = posterior_sample_2i[:,1]-posterior_sample_2i[:,3]
Deact_test_2i= posterior_sample_2i[:,2]-posterior_sample_2i[:,3]
parameter_norm = hcat(Act_test_serum,Deact_test_serum,Act_test_2i,Deact_test_2i)
writedlm("./Data/Analysis_output/posterior parameters.txt", parameter_norm, '\t')

#Simulate the distributions from the sampled posteriors
sim_dist_store_serum = Array{Float64}(250,100)
for i in 1:100
	dist_serum = generate_single_simulation_samples_m1(posterior_sample_serum[i,:],250,data_scaling_RNAseq,cell_cycle_model)
	sim_dist_store_serum[:,i] = dist_serum
end

sim_dist_store_2i = Array{Float64}(250,100)
for i in 1:100
	dist_2i = generate_single_simulation_samples_m1(posterior_sample_2i[i,:],250,data_scaling_RNAseq,cell_cycle_model)
	sim_dist_store_2i[:,i] = dist_2i
end
writedlm("./Data/Analysis_output/serum_simulations.txt", sim_dist_store_serum, '\t')
writedlm("./Data/Analysis_output/2i_simulations.txt", sim_dist_store_2i, '\t')

#Real data distributions
data_serum = readdlm("./Data/data_files/RNAseq_serum_filtereddata.txt")
data_2i = readdlm("./Data/data_files/RNAseq_2i_filtereddata.txt")
nanog_data_serum = Get_single_gene_RNAseq_data(gene_name,data_serum)
nanog_data_2i = Get_single_gene_RNAseq_data(gene_name,data_2i)
writedlm("./Data/Analysis_output/serum_truedata_nanog.txt", nanog_data_serum, '\t')
writedlm("./Data/Analysis_output/2i_truedata_nanog.txt", nanog_data_2i, '\t')
