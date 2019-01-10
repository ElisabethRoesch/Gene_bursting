#testing for basic understanding

using DifferentialEquations

base_model = @reaction_network CC begin
          c1, G0 --> G1
          c2, G1 --> G0
          c3, G1 --> G1 + P
          c4, P --> 0
end c1 c2 c3 c4

parameter_input = [0.02,.3,2.1,1.1]
tspan = (1.,20.)
          c1 = 10^parameter_input[1]
          c2 = 10^parameter_input[2]
          c3 = 1.0
          c4 = 10^parameter_input[3]
          c5 = parameter_input[4]
          #model_rates = p of DifferentialEquations
          model_rates = [c1, c2, c3, c4]
          G=2
          #Define starting state
          G1 = round(Int64,(G*c1)/(c1+c2))
          G0 = G - G1
          P = round(Int64,((c3/c4)*((G*c1)/(c1+c2))))
          u0 = [G0, G1, P]
          print(u0)
          #Define model
          prob = DiscreteProblem(u0, tspan, model_rates)
          a =solve(prob)
          jump_prob = JumpProblem(prob, Direct(), base_model)
          sol = solve(jump_prob, FunctionMap())

#no jumps
t_to_next_reaction =Array{Float64,1}(length(sol.t)-1)
t_to_next_reaction
for i in 1:length(sol.t)-1
          t_to_next_reaction[i] = sol.t[i+1]-sol.t[i]
end
string("q","q")
title=string("times to next reaction. total reactions: ",string(length(sol.t)))
plot!(reverse(sort(t_to_next_reaction)), xlab="a", ylab="b", label="Actication rate = 0.02" ,linewidth = 4)
t002=reverse(sort(t_to_next_reaction))
t02=reverse(sort(t_to_next_reaction))
t2=reverse(sort(t_to_next_reaction))

using Plots
plot(a,label =["G0","G1", "P"])
#with jumps
jump_prob_save = plot(sol, label =["G0","G1", "P"])
savefig(jump_prob_save,"jump_prob_save_2.jpeg")
#reaction_times ploted and count total reactions. run it mulitple times-> diff sols
