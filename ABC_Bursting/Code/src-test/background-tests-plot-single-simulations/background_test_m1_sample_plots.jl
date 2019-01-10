#testing for basic understanding
using DifferentialEquations
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
slow = [1.1,2.3,1.2,1.1]
parameter_input = slow
n_samples=200
c1 = 10^parameter_input[1]
c2 = 10^parameter_input[2]
c3 = 1.0
c4 = 10^parameter_input[3]
c5 = parameter_input[4]
model_rates = [c1, c2, c3, c4]
time = 4000.0
tspan= (0.0,time)
#Define starting state
G1 = round(Int64,(2*c1)/(c1+c2))
G0 = 2 - G1
P = round(Int64,((c3/c4)*((2*c1)/(c1+c2))))
u0 = [G0, G1, P]
#Define model
prob = DiscreteProblem(u0, tspan, model_rates)
jump_prob = JumpProblem(prob,Direct(),base_model)
sol = solve(jump_prob, FunctionMap())
#Simulate multiple cells and take end point
max_time = sol.t[end]
# sample_points = StatsBase.sample(0.0:max_time,n_samples)
sample_points = StatsBase.sample(0.0:max_time,n_samples)
output = sol(sample_points)
output = [data_scaling_RNAseq(i, c5) for i in output]
r=range(1,n_samples)
using Plots
b=plot(sol,xlab="time",ylab="product abundance", title="slow: sol of jumpproblem")

a=plot(r,output,title ="slow: product abundance of 100 cells/timepoints converted by scaling factor",xlab="cell/timepoint",ylab="product abundance",label =["P"])
title =string("simulated product abundance of ", n_samples," cells")
# savefig(a,"slowa.png")
# savefig(b,"slowb.png")
outputfast =output
histogram(outputfast,label=["fast product abundance"],title=title, xlab="product abundance", ylab = "frequency",nbins=20)
histogram(output,label=["product abundance"],title=title, xlab="product abundance", ylab = "frequency")




length(sol.t)
step_size = round(Int,length(sol.t)/n_samples, RoundDown)
sol.t[1:step_size:end]
