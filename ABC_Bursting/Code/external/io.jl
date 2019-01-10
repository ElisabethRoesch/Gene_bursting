import Distributions
import StatsBase
import Base: write

export
    # abstract types
    ABCInput,
    ABCOutput,

    # ABCInput structs
    ABCRejectionInput,
    ABCSMCInput,

    # ABCOutput structs
    ABCRejectionOutput,
    ABCSMCOutput,

    # Additional (mutable) structs
    ABCSMCTracker,

    # Read/write access
    read_rejection_output,
    read_smc_output

abstract type ABCInput end

struct ABCRejectionInput <: ABCInput
    n_params::Int64
    n_particles::Int64
    threshold::Float64
    priors::Vector{Distributions.ContinuousUnivariateDistribution}
    distance_function::Function
    data_generating_function::Function
    # arguments: parameter_vector,, args...
end


struct ABCSMCInput <: ABCInput
    n_params::Int64
    n_particles::Int64
    threshold_schedule::Vector{Float64}
    priors::Vector{Distributions.ContinuousUnivariateDistribution}
    distance_function::Function
    data_generating_function::Function
    # arguments: parameter_vector, args...
end


mutable struct ABCSMCTracker
    n_params::Int64
    n_accepted::Vector{Int64}
    n_tries::Vector{Int64}
    threshold_schedule::Vector{Float64}
    population::Vector{Matrix{Float64}}
    distances::Vector{Vector{Float64}}
    weights::Vector{StatsBase.Weights}
    priors::Vector{Distributions.ContinuousUnivariateDistribution}
    distance_function::Function
    data_generating_function::Function
end


abstract type ABCOutput end


struct ABCRejectionOutput <: ABCOutput
    n_params::Int64
    n_accepted::Int64
    n_tries::Int64
    threshold::Float64
    population::Matrix{Float64}  # size is n_params by n_success
    distances::Vector{Float64}
    weights::StatsBase.Weights
end


struct ABCSMCOutput <: ABCOutput
    n_params::Int64
    n_accepted::Vector{Int64}
    n_tries::Vector{Int64}
    threshold_schedule::Vector{Float64}
    population::Vector{Matrix{Float64}}
    distances::Vector{Vector{Float64}}
    weights::Vector{StatsBase.Weights}
end



function write_scalar(
        stream::IO,
        x::Real;
        separator="\n",
        )
    write(stream, string(x), separator)
end


function write_vector(
        stream::IO,
        vector::Vector{R};
        separator=" ",
        finalcharacter="\n"
        ) where {
        R<:Real,
        }
    n = length(vector)
    for i in 1:(n - 1)
        write_scalar(stream, vector[i]; separator=separator)
    end
    write_scalar(stream, vector[n]; separator=finalcharacter)
end


function write_matrix(
        stream::IO,
        matrix::Matrix{R};
        separator=" ",
        newline=true,
        ) where {
        R<:Real,
        }
    n = size(matrix, 2)
    for i in 1:(n - 1)
        write_vector(stream, matrix[:, i]; separator=separator, finalcharacter="\n")
    end
    finalcharacter = newline ? "\n" : ""
    write_vector(stream, matrix[:, n]; separator=separator, finalcharacter=finalcharacter)
end


function Base.write(stream::IO, output::ABCRejectionOutput)
    write_scalar(stream, output.n_params)
    write_scalar(stream, output.n_accepted)
    write_scalar(stream, output.n_tries)
    write_scalar(stream, output.threshold)
    write_matrix(stream, output.population)
    write_vector(stream, output.distances)
    write_vector(stream, output.weights.values)

    flush(stream)
end


function Base.write(stream::IO, output::ABCSMCOutput)
    write_scalar(stream, output.n_params)
    write_vector(stream, output.n_accepted)
    write_vector(stream, output.n_tries)
    write_vector(stream, output.threshold_schedule)
    for i in 1:length(output.threshold_schedule)
        write(stream, "-----\n")
        write_scalar(stream, output.n_accepted[i])
        write_scalar(stream, output.n_tries[i])
        write_scalar(stream, output.threshold_schedule[i])
        write_matrix(stream, output.population[i])
        write_vector(stream, output.distances[i])
        write_vector(stream, output.weights[i].values)
    end

    flush(stream)
end


function write(stream::IO, output::ABCOutput)
    return Base.write(stream, output)
end


function read_rejection_output(filepath::AbstractString)
    input_file = open(filepath, "r")

    try
        n_params = parse(Int64, readline(input_file))
        n_accepted = parse(Int64, readline(input_file))
        n_tries = parse(Int64, readline(input_file))
        threshold = parse(Float64, readline(input_file))

        population = zeros(n_params, n_accepted)
        for i in 1:n_accepted
            population[:, i] = parse.(Float64, split(readline(input_file)))
        end

        distances = parse.(Float64, split(readline(input_file)))

        wts = parse.(Float64, split(readline(input_file)))
        weights = StatsBase.Weights(wts)

        return ABCRejectionOutput(n_params,
                                  n_accepted,
                                  n_tries,
                                  threshold,
                                  population,
                                  distances,
                                  weights,
                                  )
    finally
        close(input_file)
    end
end


function read_smc_output(filepath::AbstractString)
    input_file = open(filepath, "r")

    try
        n_params = parse(Int64, readline(input_file))
        n_accepted = parse.(Int64, split(readline(input_file)))
        n_tries = parse.(Int64, split(readline(input_file)))
        threshold_schedule = parse.(Float64, split(readline(input_file)))

        @assert length(n_accepted) == length(n_tries) == length(threshold_schedule)

        population = Matrix{Float64}[]
        distances = Vector{Float64}[]
        weights = StatsBase.Weights[]

        for i in 1:length(threshold_schedule)
            separator = readline(input_file)

            @assert parse(Int64, readline(input_file)) == n_accepted[i]
            @assert parse(Int64, readline(input_file)) == n_tries[i]
            @assert parse(Float64, readline(input_file)) == threshold_schedule[i]

            push!(population, zeros(n_params, n_accepted[i]))
            for j in 1:n_accepted[i]
                particle = parse.(Float64, split(readline(input_file)))
                population[i][:, j] = particle
            end

            push!(distances, parse.(Float64, split(readline(input_file))))

            wts = parse.(Float64, split(readline(input_file)))
            push!(weights, StatsBase.Weights(wts))
        end

        return ABCSMCOutput(n_params,
                            n_accepted,
                            n_tries,
                            threshold_schedule,
                            population,
                            distances,
                            weights,
                            )
    finally
        close(input_file)
    end
end
