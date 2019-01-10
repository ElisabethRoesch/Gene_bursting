module ABaCus
import Distributions
import StatsBase


export generate_kernels, generate_parameters,
       ABCrejection,
       initialiseABCSMC, iterateABCSMC!, ABCSMC

include("io.jl")


function generate_kernels(
        population::Matrix{F},
        priors::Vector{D},
        ) where {
        F<:AbstractFloat,
        D<:Distributions.ContinuousUnivariateDistribution,
        }
    n_params = size(population, 1)
    n_particles = size(population, 2)

    stds = std(population, 2)[:]
    lowers = minimum.(priors)
    uppers = maximum.(priors)

    CUD = Distributions.ContinuousUnivariateDistribution
    kernels = Matrix{CUD}(n_params, n_particles)
    for j in 1:n_particles
        means = population[:, j]
        for i in 1:n_params
            kernels[i, j] = Distributions.TruncatedNormal(means[i],
                                                          stds[i]*sqrt(2),
                                                          lowers[i],
                                                          uppers[i],
                                                          )
        end
    end

    return kernels
end


"""
A function which generates a parameter vector from the prior distribution.
"""
function generate_parameters(
        priors::Vector{D},
        ) where {
        D<:Distributions.ContinuousUnivariateDistribution
        }
    n_dims = length(priors)
    parameters = zeros(n_dims)

    weight = 1.
    for i in 1:n_dims
        @inbounds parameters[i] = rand(priors[i])
        weight *= Distributions.pdf(priors[i], parameters[i])
    end

    return parameters, weight
end


function generate_parameters(
        priors::Vector{D1},
        old_parameters::Matrix{F},
        old_weights::StatsBase.Weights,
        kernels::Matrix{D2},
        ) where {
        D1, D2<:Distributions.ContinuousUnivariateDistribution,
        F<:AbstractFloat,
        }
    n_params = length(priors)

    # ADD DimensionMismatch THROWS SO @inbounds CAN BE USED?

    # the kernels must be centered around the old particles
    # and truncated to the priors.

    particle = StatsBase.sample(indices(old_parameters, 2), old_weights)
    perturbed_parameters = rand.(kernels[:, particle])

    numerator = 1.0
    for i in 1:n_params
        numerator *= Distributions.pdf(priors[i], perturbed_parameters[i])
    end

    denominator = 0.0
    for k in eachindex(old_weights)
        # calculate the total kernel
        kernel = 1.0
        for j in 1:n_params
            kernel *= Distributions.pdf(kernels[j, k], perturbed_parameters[j])
        end
        denominator += old_weights[k] * kernel

        # weight normalisation---for numerical stability if nothing else
        denominator *= old_weights.sum
    end

    weight = numerator / denominator

    return perturbed_parameters, weight
end


function normalise(
        weights::StatsBase.AbstractWeights;
        tosum = 1.0,
        )
    WeightType = typeof(weights)
    weights = WeightType(weights.values .* (tosum / sum(weights.values)), tosum)

    return weights
end


"""
The ABC rejection sampler.
"""
function ABCrejection(
        input::ABCRejectionInput,
        reference_data,
        simulation_args...;
        out_stream::IO = STDOUT,
        )
    # initialise
    n_tries = 0
    n_accepted = 0
    accepted_parameters = zeros(input.n_params, input.n_particles)
    accepted_distances = zeros(input.n_particles)
    weights = ones(input.n_particles)
    count = n_tries+1

    # simulate
    while n_accepted < input.n_particles
        parameters, weight = generate_parameters(input.priors)
        simulated_data = input.data_generating_function(parameters,
                                                        simulation_args...
                                                        )
        distance = input.distance_function(reference_data, simulated_data)
        n_tries += 1

        if distance < input.threshold
            n_accepted += 1
            accepted_parameters[:, n_accepted] = parameters
            accepted_distances[n_accepted] = distance
            weights[n_accepted] = weight
        end

        if n_tries == count
            println("Number of trials = $n_tries")
            println("Number accepted = $n_accepted")
            count = count+1000
        end


    end

    weights = weights ./ sum(weights)

    # output
    output = ABCRejectionOutput(input.n_params,
                                n_accepted,
                                n_tries,
                                input.threshold,
                                accepted_parameters,
                                accepted_distances,
                                StatsBase.Weights(weights, 1.0),
                                )

    write(out_stream, output)

    return output
end


"""
Simulate the first round of ABC-SMC.
Output: a ABCSMCTracker.
"""
function initialiseABCSMC(
        input::ABCSMCInput,
        reference_data,
        simulation_args...
        )
    # the first run is an ABC rejection simulation
    rejection_input = ABCRejectionInput(input.n_params,
                                        input.n_particles,
                                        input.threshold_schedule[1],
                                        input.priors,
                                        input.distance_function,
                                        input.data_generating_function,
                                        )
    rejection_output = ABCrejection(rejection_input,
                                    reference_data,
                                    simulation_args...
                                    )

    tracker =  ABCSMCTracker(input.n_params,
                             [rejection_output.n_accepted],
                             [rejection_output.n_tries],
                             [rejection_output.threshold],
                             [rejection_output.population],
                             [rejection_output.distances],
                             [rejection_output.weights],
                             input.priors,
                             input.distance_function,
                             input.data_generating_function,
                             )

    return tracker
end


"""
Simulate one round of ABC-SMC (not the first)
"""
function iterateABCSMC!(
        tracker::ABCSMCTracker,
        threshold::AbstractFloat,
        n_toaccept::Integer,
        reference_data,
        simulation_args...
        )
    # initialise
    push!(tracker.n_accepted, 0)
    push!(tracker.n_tries, 0)
    if threshold > tracker.threshold_schedule[end]
        println("Warning: current threshold less strict than previous one.")
    end
    push!(tracker.threshold_schedule, threshold)
    push!(tracker.population, zeros(tracker.population[end]))
    push!(tracker.distances, zeros(tracker.distances[end]))
    push!(tracker.weights,
          StatsBase.Weights(ones(tracker.weights[end].values)))

    kernels = generate_kernels(tracker.population[end-1], tracker.priors)

    # simulate
    while tracker.n_accepted[end] < n_toaccept
        parameters, weight = generate_parameters(tracker.priors,
                                                 tracker.population[end-1],
                                                 tracker.weights[end-1],
                                                 kernels,
                                                 )
        simulated_data = tracker.data_generating_function(parameters,
                                                          simulation_args...
                                                          )
        distance = tracker.distance_function(reference_data, simulated_data)
        tracker.n_tries[end] += 1

        if distance < threshold
            tracker.n_accepted[end] += 1
            n_accepted = tracker.n_accepted[end]
            tracker.population[end][:, n_accepted] = parameters
            tracker.distances[end][n_accepted] = distance
            tracker.weights[end].values[n_accepted] = weight
        end
    end

    tracker.weights[end] = deepcopy(normalise(tracker.weights[end], tosum=1.0))

    return tracker
end


"""
Simulate an entire ABC-SMC run.
"""
function ABCSMC(
        input::ABCSMCInput,
        reference_data,
        simulation_args...;
        out_stream::IO = STDOUT,
        )
    n_toaccept = input.n_particles

    tracker = initialiseABCSMC(input, reference_data, simulation_args...)

    for i in 2:length(input.threshold_schedule)
        threshold = input.threshold_schedule[i]
        iterateABCSMC!(tracker,
                       threshold,
                       input.n_particles,
                       reference_data,
                       simulation_args...
                       )
        output = ABCSMCOutput(input.n_params,
                              tracker.n_accepted,
                              tracker.n_tries,
                              tracker.threshold_schedule,
                              tracker.population,
                              tracker.distances,
                              tracker.weights,
                              )

        write(out_stream, output)
    end

    return ABCSMCOutput(input.n_params,
                        tracker.n_accepted,
                        tracker.n_tries,
                        tracker.threshold_schedule,
                        tracker.population,
                        tracker.distances,
                        tracker.weights,
                        )
end


end # module
