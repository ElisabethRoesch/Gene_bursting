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

slow = [1.,1.,-1.,1.]
slow = [1.1,0.7,-1.,1.]

parameter_input = slow
n_samples=250
cycle_time_vec = [2000.,2000.]

#act,de-act,deg
max_time = 40000.0
sample_points = StatsBase.sample(0.0:max_time,n_samples)
time_1,time_2 = cycle_time_vec  # time1/2: time it takes for first/second part of cell cycle
c1 = 10^parameter_input[1]
c2 = 10^parameter_input[2]
c3 = 1.0
c4 = 10^parameter_input[3]
c5 = parameter_input[4]
p = [c1,c2,c3,c4]
#Define initial starting state (starting in 2-state system)
G = 2
G1 = round(Int64,(2*c1)/(c1+c2))
G0 = G - G1
P = round(Int64,((c3/c4)*((2*c1)/(c1+c2))))
starting_state = [G0,G1,P]
time_range = (0.0,time_1)

output_samples = []
start_time = 0.
end_time = 0.
sample_points = StatsBase.sample(0.0:max_time,n_samples)
sols= []
while end_time < max_time
            #Cell cycle part 1
            prob = DiscreteProblem(starting_state,time_range,p)
            jump_prob = JumpProblem(prob,Direct(),base_model)
            sol = solve(jump_prob,FunctionMap())
            push!(sols,sol)
            start_time = sol.t[1]
            end_time = sol.t[end]
            sub_sample_points = sample_points[[x >= start_time for x in sample_points]]
            sub_sample_points = sub_sample_points[[x < end_time for x in sub_sample_points]]
            if length(sub_sample_points) >= 1
                output = sol(sub_sample_points)
                print(output)
                output = [data_scaling_RNAseq(i, c5) for i in output]
                output_samples = vcat(output_samples,output)
            end
            #Cell cycle part 2
            starting_state = sol.u[end] #Redefine the starting state based on end point of previous simulation
            starting_state[1] = starting_state[1]*2 #Double number of genes, keeping their current status of activitity (on/off)
            starting_state[2] = starting_state[2]*2
            time_range= (end_time,end_time+time_2)
            prob = DiscreteProblem(starting_state,time_range,p)
            jump_prob = JumpProblem(prob,Direct(),base_model)
            sol = solve(jump_prob,FunctionMap())
            start_time = sol.t[1]
            end_time = sol.t[end]
            sub_sample_points = sample_points[[x >= start_time for x in sample_points]]
            sub_sample_points = sub_sample_points[[x < end_time for x in sub_sample_points]]
            if length(sub_sample_points) >= 1
                output = sol(sub_sample_points)
                output = [data_scaling_RNAseq(i, c5) for i in output]
                output_samples = vcat(output_samples,output)
            end
            ##Cell cycle part 1 starting state redefined
            starting_state = sol.u[end] #Redefine the starting state based on end point of previous simulation
            if starting_state[1] <= 1
                starting_state[1] = 0
                starting_state[2] = 2
            elseif starting_state[2] <= 1
                starting_state[1] = 2
                starting_state[2] = 0
            else
                starting_state[1] = 1
                starting_state[2] = 1
            end
            push!(sols,sol)
            time_range=(end_time,end_time+time_1)
end

output = float([i for i in output_samples])

r=range(1,n_samples)
using Plots
# b=scatter(sols,xlab="time",ylab="product abundance", title="slow: sol of jumpproblem")

a=plot(r,output,title ="Product abundance of 200 cells/timepoints converted by scaling factor",xlab="cell/timepoint",ylab="product abundance",label =["P"])
title =string("simulated product abundance of ", n_samples," cells")
# savefig(a,"slowa.png")
# savefig(b,"slowb.png")

histogram(output,label=["fast product abundance"],title=title, xlab="product abundance", ylab = "frequency",nbins=20)
histogram(output,label=["product abundance"],title=title, xlab="product abundance", ylab = "frequency")




t_to_next_reaction1 =Array{Float64,1}(length(sols[1].t)-1)
t_to_next_reaction2 =Array{Float64,1}(length(sols[2].t)-1)
t_to_next_reaction3 =Array{Float64,1}(length(sols[3].t)-1)
t_to_next_reaction4 =Array{Float64,1}(length(sols[4].t)-1)

for i in 1:length(sols[1].t)-1
          t_to_next_reaction1[i] = sols[1].t[i+1]-sols[1].t[i]
end
for i in 1:length(sols[2].t)-1
          t_to_next_reaction2[i] = sols[2].t[i+1]-sols[2].t[i]
end
for i in 1:length(sols[3].t)-1
          t_to_next_reaction3[i] = sols[3].t[i+1]-sols[3].t[i]
end
for i in 1:length(sols[4].t)-1
          t_to_next_reaction4[i] = sols[4].t[i+1]-sols[4].t[i]
end
string(string(length(sols[1].t)),",",1)
title=string("times to next reaction. total reactions: ",length(sols[1].t), ", ",length(sols[2].t), ", ",length(sols[3].t),", ",length(sols[4].t) )
plot(reverse(sort(t_to_next_reaction1)), title=title,label="CC1 part 1")
plot!(reverse(sort(t_to_next_reaction2)),label="CC1 part 2")
plot!(reverse(sort(t_to_next_reaction3)),label="CC2 part 1")
plot!(reverse(sort(t_to_next_reaction4)),label="CC2 part 2")

jump_prob_save = plot(sols[1], label =["G0","G1", "P"])
jump_prob_save = plot(sols[2], label =["G0","G1", "P"])
jump_prob_save = plot(sols[3], label =["G0","G1", "P"])
jump_prob_save = plot(sols[4], label =["G0","G1", "P"])
