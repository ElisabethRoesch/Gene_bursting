############################################# file_description ###############################################
#this file contains     1. the base_model as a reaction_network of DifferentialEquations.
#                       2. the function generate_single_simulation_m1:
                                    # in:  sampled parameter, tspan
                                    # out: solution of solve of DifferentialEquations
#                       3. the function generate_single_simulation_samples_m1:
                                    # in:  sampled parameter, #timepoints, data_scaling_function
                                    # out: production estimation
##############################################################################################################

############################################### base_model ###################################################
base_model = @reaction_network CC begin
c1, G0 --> G1
c2, G1 --> G0
c3, G1 --> G1 + P
c4, P --> 0
end c1 c2 c3 c4
##############################################################################################################


################################################ m_1 ##########################################################
#Simulated product throughout time
function generate_single_simulation_m1(
            parameter_input::Vector{Float64},
            tspan::Tuple{Float64,Float64}
            )
            #she traslates the parameters sampled from the prior into folds
            c1 = 10^parameter_input[1]
            c2 = 10^parameter_input[2]
            c3 = 1.0
            c4 = 10^parameter_input[3]
            c5 = 1.
            model_rates = [c1, c2, c3, c4]
            #Define starting state
            G1 = round(Int64,(2*c1)/(c1+c2))
            G0 = 2 - G1
            P = round(Int64,((c3/c4)*((2*c1)/(c1+c2))))
            u0 = [G0, G1, P]
            #Define model
            prob = DiscreteProblem(u0, tspan, model_rates)
            jump_prob = JumpProblem(prob, Direct(), base_model)
            sol = solve(jump_prob, FunctionMap())
    return sol
end
#Simulated product throughout time and randomly samples to produce a distribution
function generate_single_simulation_samples_m1(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )
            c1 = 10^parameter_input[1]
            c2 = 10^parameter_input[2]
            c3 = 1.0
            c4 = 10^parameter_input[3]
            c5 = 1.
            model_rates = [c1, c2, c3, c4]
            time = 40000.0
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
            #change to grid?
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            output = sol(sample_points)
            #dimesion reducting to product (3rd species)
            output = [data_scaling_function(i, c5) for i in output]
    return output
end
function generate_single_simulation_samples_m1_25000(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )
            c1 = 10^parameter_input[1]
            c2 = 10^parameter_input[2]
            c3 = 1.0
            c4 = 10^parameter_input[3]
            c5 = 1.0
            model_rates = [c1, c2, c3, c4]
            time = 25000.0
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
            #change to grid?
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            output = sol(sample_points)
            #dimesion reducting to product (3rd species)
            output = [data_scaling_function(i, c5) for i in output]
    return output
end
function generate_single_simulation_samples_m1_50000(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )
            c1 = 10^parameter_input[1]
            c2 = 10^parameter_input[2]
            c3 = 1.0
            c4 = 10^parameter_input[3]
            c5 = 1.0
            model_rates = [c1, c2, c3, c4]
            time = 50000.0
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
            #change to grid?
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            output = sol(sample_points)
            #dimesion reducting to product (3rd species)
            output = [data_scaling_function(i, c5) for i in output]
    return output
end
#Simulated product throughout time and randomly samples to produce a distribution
function generate_single_simulation_samples_m1_100000(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )
            c1 = 10^parameter_input[1]
            c2 = 10^parameter_input[2]
            c3 = 1.0
            c4 = 10^parameter_input[3]
            c5 = 1.
            model_rates = [c1, c2, c3, c4]
            time = 100000.0
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
            #change to grid?
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            output = sol(sample_points)
            #dimesion reducting to product (3rd species)
            output = [data_scaling_function(i, c5) for i in output]
    return output
end
function generate_single_simulation_samples_m1_125000(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )
            c1 = 10^parameter_input[1]
            c2 = 10^parameter_input[2]
            c3 = 1.0
            c4 = 10^parameter_input[3]
            c5 = 1.
            model_rates = [c1, c2, c3, c4]
            time = 125000.0
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
            #change to grid?
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            output = sol(sample_points)
            #dimesion reducting to product (3rd species)
            output = [data_scaling_function(i, c5) for i in output]
    return output
end

#Simulated product throughout time and randomly samples to produce a distribution
function generate_single_simulation_samples_m1_200000(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )
            c1 = 10^parameter_input[1]
            c2 = 10^parameter_input[2]
            c3 = 1.0
            c4 = 10^parameter_input[3]
            c5 = 1.
            model_rates = [c1, c2, c3, c4]
            time = 200000.0
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
            #change to grid?
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            output = sol(sample_points)
            #dimesion reducting to product (3rd species)
            output = [data_scaling_function(i, c5) for i in output]
    return output
end
function generate_single_simulation_samples_m1_250000(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )
            c1 = 10^parameter_input[1]
            c2 = 10^parameter_input[2]
            c3 = 1.0
            c4 = 10^parameter_input[3]
            c5 = 1.
            model_rates = [c1, c2, c3, c4]
            time = 250000.0
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
            #change to grid?
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            output = sol(sample_points)
            #dimesion reducting to product (3rd species)
            output = [data_scaling_function(i, c5) for i in output]
    return output
end
function generate_single_simulation_samples_m1_500000(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )
            c1 = 10^parameter_input[1]
            c2 = 10^parameter_input[2]
            c3 = 1.0
            c4 = 10^parameter_input[3]
            c5 = 1.0
            model_rates = [c1, c2, c3, c4]
            time = 500000.0
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
            #change to grid?
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            output = sol(sample_points)
            #dimesion reducting to product (3rd species)
            output = [data_scaling_function(i, c5) for i in output]
    return output
end
###################################################################################################


#########################################m_8#######################################################


function generate_single_simulation_m8(
            parameter_input::Vector{Float64},
            time_range::Tuple{Float64,Float64},
            )
            #Define starting state
           c1 = 10^parameter_input[1]  #Activation
           c2 = 10^parameter_input[2]  #Deactivation
           c3 = 1.0                    #Expression
           c4 = 10^parameter_input[3]  #Degradation
           k = 10^parameter_input[4]     #Feedback
           c6 = 1.   #Scale factor
           #Define starting state at stationary distribution
           G = 2.0
           P = (sqrt((4*c1*c2*c3*c4*G*k)+(c1*c4*k+c2*c4*k)^2)-c1*c4*k-c2*c4*k)/(2*c2*c4)
           G1 = round(Int64,(
                                                (c1*G*k)/((c1*k)+(c2*(k+P)))
           ))
           P = round(Int64,(P))
           if P < 0
                                                P = 0
           else
                                                P = P
           end
           if G1 < 0
                                                G1 = 0
           else
                                                G1 = G1
           end
           G0 = 2 - G1
           starting_state = [G0,G1,P]
           time = 40000.0
           time_range= (0.0,time)
           #Define model
           #Set rate parameters
                     rate1 = function (u,p,t)
                                return (c1)*u[1]*(1/(1+(u[3]/k)))
                     end

                     affect! = function (integrator)
                                integrator.u[1] -= 1
                                integrator.u[2] += 1
                     end
                     jump1 = DifferentialEquations.ConstantRateJump(rate1,affect!)

                     rate2 = function (u,p,t)
                                return (c2)*u[2]
                     end

                     affect! = function (integrator)
                                integrator.u[1] += 1
                                integrator.u[2] -= 1
                     end
                     jump2 = DifferentialEquations.ConstantRateJump(rate2,affect!)

                     rate3 = function (u,p,t)
                                return (c3)*u[2]
                     end
                     affect! = function (integrator)
                                integrator.u[2] += 0
                                integrator.u[3] += 1
                     end
                     jump3 = DifferentialEquations.ConstantRateJump(rate3,affect!)

                     rate4 = function (u,p,t)
                                return (c4)*u[3]
                     end

                     affect! = function (integrator)
                                integrator.u[3] -= 1
                     end
                     jump4 = DifferentialEquations.ConstantRateJump(rate4,affect!)

                     prob = DifferentialEquations.DiscreteProblem(starting_state,time_range)
                     jump_prob = DifferentialEquations.JumpProblem(prob,DifferentialEquations.Direct(),jump1,jump2,jump3,jump4)

                     #Simulate multiple cells and take end point
                     distribution_data = Array{Float64}(n_samples)
                     sol = DifferentialEquations.solve(jump_prob,FunctionMap())
       return sol
end


function generate_single_simulation_samples_m8(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )
	#Define starting state
	c1 = 10^parameter_input[1]  #Activation
	c2 = 10^parameter_input[2]  #Deactivation
	c3 = 1.0                    #Expression
	c4 = 10^parameter_input[3]  #Degradation
	k = 10^parameter_input[4]     #Feedback
	c6 = 1.     #Scale factor
	#Define starting state at stationary distribution
	G = 2.0
	P = (sqrt((4*c1*c2*c3*c4*G*k)+(c1*c4*k+c2*c4*k)^2)-c1*c4*k-c2*c4*k)/(2*c2*c4)
	G1 = round(Int64,(
		(c1*G*k)/((c1*k)+(c2*(k+P)))
	))
	P = round(Int64,(P))
	if P < 0
		P = 0
	else
		P = P
	end
	if G1 < 0
		G1 = 0
	else
		G1 = G1
	end
	G0 = 2 - G1
	starting_state = [G0,G1,P]
	time = 250000.0
	time_range= (0.0,time)
	#Define model
	#Set rate parameters


                                    rate1 = function (u,p,t)
                                               return (c1)*u[1]*(1/(1+(u[3]/k)))
                                    end

                                    affect! = function (integrator)
                                              integrator.u[1] -= 1
                                              integrator.u[2] += 1
                                    end
                                    jump1 = DifferentialEquations.ConstantRateJump(rate1,affect!)

                                    rate2 = function (u,p,t)
                                              return (c2)*u[2]
                                    end

                                    affect! = function (integrator)
                                              integrator.u[1] += 1
                                              integrator.u[2] -= 1
                                    end
                                    jump2 = DifferentialEquations.ConstantRateJump(rate2,affect!)

                                    rate3 = function (u,p,t)
                                              return (c3)*u[2]
                                    end
                                    affect! = function (integrator)
                                              integrator.u[2] += 0
                                              integrator.u[3] += 1
                                    end
                                    jump3 = DifferentialEquations.ConstantRateJump(rate3,affect!)

                                    rate4 = function (u,p,t)
                                              return (c4)*u[3]
                                    end

                                    affect! = function (integrator)
                                              integrator.u[3] -= 1
                                    end
                                    jump4 = DifferentialEquations.ConstantRateJump(rate4,affect!)

                                    prob = DifferentialEquations.DiscreteProblem(starting_state,time_range)
                                    jump_prob = DifferentialEquations.JumpProblem(prob,DifferentialEquations.Direct(),jump1,jump2,jump3,jump4)

                                    #Simulate multiple cells and take end point
                                    distribution_data = Array{Float64}(n_samples)
                                    sol = DifferentialEquations.solve(jump_prob,FunctionMap())




	max_time = sol.t[end]
	for i in 1:n_samples
		outcome = round.(Int64,(sol(i*(max_time/(n_samples+1)))))
		distribution_data[i,1] = data_scaling_function(outcome, c6)
	end
	return distribution_data
end

#################################################################################################################################################

########################################################## m_9 #############################################################################

function generate_single_simulation_m9(
            parameter_input::Vector{Float64},
            time_range::Tuple{Float64,Float64},
            )
	#Define starting state
	c1 = 10^parameter_input[1]  #Activation
	c2 = 10^parameter_input[2]  #Deactivation
	c3 = 1.0                    #Expression
	c4 = 10^parameter_input[3]  #Degradation
	c5 = 10^parameter_input[4]     #Feedback baseline
	k = 10^parameter_input[5]     #Feedback k
	c6 = 1.     #Scale factor


	#Define starting state at stationary distribution
	G = 2

	P =(
	        (sqrt(4*c1*c3*c4*G*k*(c1+c2+c5)+(c1*c4*k-c1*c3*G+c2*c4*k)^2)
	          +
	         (c1*c3*G-c1*c4*k-c2*c4*k)
	         )
	           /
	(2*c4*(c1+c2+c5))
	)

	G1 = round(Int64,(
	(c1*G*(k+P))
	/
	(c1*(k+P)+c2*(k+P)+c5*P)
	))
	P = round(Int64,(P))
	if P < 0
		P = 0
	else
		P = P
	end
	if G1 < 0
		G1 = 0
	else
		G1 = G1
	end
	G0 = 2 - G1
	starting_state = [G0,G1,P]
	#Set rate parameters
	rate = (t,u) -> (c1)*u[1]
	affect! = function (integrator)
	integrator.u[1] -= 1
	integrator.u[2] += 1
	end
	jump1 = DifferentialEquations.ConstantRateJump(rate,affect!)
	rate = (t,u) -> (c2*u[2])+(u[2]*c5*((u[3]/k)/(1+(u[3]/k))))
	affect! = function (integrator)
	integrator.u[1] += 1
	integrator.u[2] -= 1
	end
	jump2 = DifferentialEquations.ConstantRateJump(rate,affect!)
	rate = (t,u) -> (c3)*u[2]
	affect! = function (integrator)
	integrator.u[2] += 0
	integrator.u[3] += 1
	end
	jump3 = DifferentialEquations.ConstantRateJump(rate,affect!)
	rate = (t,u) -> (c4)*u[3]
	affect! = function (integrator)
	integrator.u[3] -= 1
	end
	jump4 = DifferentialEquations.ConstantRateJump(rate,affect!)
	prob = DifferentialEquations.DiscreteProblem(starting_state,time_range)
	jump_prob = DifferentialEquations.JumpProblem(prob,DifferentialEquations.Direct(),jump1,jump2,jump3,jump4)
	sol = DifferentialEquations.solve(jump_prob,DifferentialEquations.Discrete(apply_map=false))
	return sol
end


function generate_single_simulation_samples_m9(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )

            #Define starting state
            c1 = 10^parameter_input[1]  #Activation
            c2 = 10^parameter_input[2]  #Deactivation
            c3 = 1.0                    #Expression
            c4 = 10^parameter_input[3]  #Degradation
            c5 = 10^parameter_input[4]     #Feedback baseline
            k =  10^parameter_input[5]     #Feedback k
            c6 = 0.977   #Scale factor


            #Define starting state at stationary distribution

            G = 2

            P =(
                      (
                                sqrt(4*c1*c3*c4*G*k*(c1+c2+c5)+(c1*c4*k-c1*c3*G+c2*c4*k)^2)
                                +
                                (c1*c3*G-c1*c4*k-c2*c4*k)
                      )
                      /
                      (2*c4*(c1+c2+c5))
            )

            G1 = round(Int64,(

            (c1*G*(k+P))
            /
            (c1*(k+P)+c2*(k+P)+c5*P)
            ))

            P = round(Int64,(P))
            if P < 0
                      P = 0
            else
                      P = P
            end

            if G1 < 0
                      G1 = 0
            else
                      G1 = G1
            end
            G0 = 2 - G1
            starting_state = [G0,G1,P]
            time = 250000.0
            time_range= (0.0,time)
            #Define model
            #Set rate parameters
            rate1 = function (u,p,t)
                                                return (c1)*u[1]
            end

            affect! = function (integrator)
                      integrator.u[1] -= 1
                      integrator.u[2] += 1
            end
            jump1 = DifferentialEquations.ConstantRateJump(rate1,affect!)

            rate2 = function (u,p,t)
                      return (c2*u[2])+(u[2]*c5*((u[3]/k)/(1+(u[3]/k))))
            end

            affect! = function (integrator)
                      integrator.u[1] += 1
                      integrator.u[2] -= 1
            end
            jump2 = DifferentialEquations.ConstantRateJump(rate2,affect!)

            rate3 = function (u,p,t)
                      return (c3)*u[2]
             end
            affect! = function (integrator)
                      integrator.u[2] += 0
                      integrator.u[3] += 1
            end
            jump3 = DifferentialEquations.ConstantRateJump(rate3,affect!)

            rate4 = function (u,p,t)
                      return (c4)*u[3]
            end
            affect! = function (integrator)
                      integrator.u[3] -= 1
            end
            jump4 = DifferentialEquations.ConstantRateJump(rate4,affect!)

            prob = DifferentialEquations.DiscreteProblem(starting_state,time_range)
            jump_prob = DifferentialEquations.JumpProblem(prob,DifferentialEquations.Direct(),jump1,jump2,jump3,jump4)

            #Simulate multiple cells and take end point
            distribution_data = Array{Float64}(n_samples)
            sol = DifferentialEquations.solve(jump_prob,DifferentialEquations.FunctionMap())
            max_time = sol.t[end]
            for i in 1:n_samples
                      outcome = round.(Int64,(sol(i*(max_time/(n_samples+1)))))
                      distribution_data[i,1] = data_scaling_RNAseq(outcome, c6)
            end
      return distribution_data
end


###try:

function generate_single_simulation_samples_m9_no_deact(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )

            c1 = 10^parameter_input[1]  #Activation
            c2 = 0.  #Deactivation
            c3 = 1.0                    #Expression
            c4 = 10^parameter_input[2]  #Degradation
            c5 = 10^parameter_input[3]     #Feedback baseline
            k =  10^parameter_input[4]     #Feedback k
            c6 = 0.977   #Scale factor
            G = 2
            P =(
                      (
                                sqrt(4*c1*c3*c4*G*k*(c1+c2+c5)+(c1*c4*k-c1*c3*G+c2*c4*k)^2)
                                +
                                (c1*c3*G-c1*c4*k-c2*c4*k)
                      )
                      /
                      (2*c4*(c1+c2+c5))
            )

            G1 = round(Int64,(

            (c1*G*(k+P))
            /
            (c1*(k+P)+c2*(k+P)+c5*P)
            ))

            P = round(Int64,(P))
            if P < 0
                      P = 0
            else
                      P = P
            end

            if G1 < 0
                      G1 = 0
            else
                      G1 = G1
            end
            G0 = 2 - G1
            starting_state = [G0,G1,P]
            time = 250000.0
            time_range= (0.0,time)
            #Define model
            #Set rate parameters
            rate1 = function (u,p,t)
                                                return (c1)*u[1]
            end

            affect! = function (integrator)
                      integrator.u[1] -= 1
                      integrator.u[2] += 1
            end
            jump1 = DifferentialEquations.ConstantRateJump(rate1,affect!)

            rate2 = function (u,p,t)
                      return (c2*u[2])+(u[2]*c5*((u[3]/k)/(1+(u[3]/k))))
            end

            affect! = function (integrator)
                      integrator.u[1] += 1
                      integrator.u[2] -= 1
            end
            jump2 = DifferentialEquations.ConstantRateJump(rate2,affect!)

            rate3 = function (u,p,t)
                      return (c3)*u[2]
             end
            affect! = function (integrator)
                      integrator.u[2] += 0
                      integrator.u[3] += 1
            end
            jump3 = DifferentialEquations.ConstantRateJump(rate3,affect!)

            rate4 = function (u,p,t)
                      return (c4)*u[3]
            end
            affect! = function (integrator)
                      integrator.u[3] -= 1
            end
            jump4 = DifferentialEquations.ConstantRateJump(rate4,affect!)

            prob = DifferentialEquations.DiscreteProblem(starting_state,time_range)
            jump_prob = DifferentialEquations.JumpProblem(prob,DifferentialEquations.Direct(),jump1,jump2,jump3,jump4)

            #Simulate multiple cells and take end point
            distribution_data = Array{Float64}(n_samples)
            sol = DifferentialEquations.solve(jump_prob,DifferentialEquations.FunctionMap())
            max_time = sol.t[end]
            for i in 1:n_samples
                      outcome = round.(Int64,(sol(i*(max_time/(n_samples+1)))))
                      distribution_data[i,1] = data_scaling_RNAseq(outcome, c6)
            end
      return distribution_data
end



#############################################################################################################################

function generate_single_simulation_samples_m9_no_deact_half(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )

            #Define starting state
            c1 = 10^parameter_input[1]  #Activation
            c2 = 0.  #Deactivation
            c3 = 1.0                    #Expression
            c4 = 10^parameter_input[2]  #Degradation
            c5 = 10^parameter_input[3]     #Feedback baseline
            k =  10^parameter_input[4]     #Feedback k
            c6 = 0.977   #Scale factor


            #Define starting state at stationary distribution

            G = 2

            P =(
                      (
                                sqrt(4*c1*c3*c4*G*k*(c1+c2+c5)+(c1*c4*k-c1*c3*G+c2*c4*k)^2)
                                +
                                (c1*c3*G-c1*c4*k-c2*c4*k)
                      )
                      /
                      (2*c4*(c1+c2+c5))
            )

            G1 = round(Int64,(

            (c1*G*(k+P))
            /
            (c1*(k+P)+c2*(k+P)+c5*P)
            ))

            P = round(Int64,(P))
            if P < 0
                      P = 0
            else
                      P = P
            end

            if G1 < 0
                      G1 = 0
            else
                      G1 = G1
            end
            G0 = 2 - G1
            starting_state = [G0,G1,P]
            time = 125000.0
            time_range= (0.0,time)
            #Define model
            #Set rate parameters
            rate1 = function (u,p,t)
                                                return (c1)*u[1]
            end

            affect! = function (integrator)
                      integrator.u[1] -= 1
                      integrator.u[2] += 1
            end
            jump1 = DifferentialEquations.ConstantRateJump(rate1,affect!)

            rate2 = function (u,p,t)
                      return (c2*u[2])+(u[2]*c5*((u[3]/k)/(1+(u[3]/k))))
            end

            affect! = function (integrator)
                      integrator.u[1] += 1
                      integrator.u[2] -= 1
            end
            jump2 = DifferentialEquations.ConstantRateJump(rate2,affect!)

            rate3 = function (u,p,t)
                      return (c3)*u[2]
             end
            affect! = function (integrator)
                      integrator.u[2] += 0
                      integrator.u[3] += 1
            end
            jump3 = DifferentialEquations.ConstantRateJump(rate3,affect!)

            rate4 = function (u,p,t)
                      return (c4)*u[3]
            end
            affect! = function (integrator)
                      integrator.u[3] -= 1
            end
            jump4 = DifferentialEquations.ConstantRateJump(rate4,affect!)

            prob = DifferentialEquations.DiscreteProblem(starting_state,time_range)
            jump_prob = DifferentialEquations.JumpProblem(prob,DifferentialEquations.Direct(),jump1,jump2,jump3,jump4)

            #Simulate multiple cells and take end point
            distribution_data = Array{Float64}(n_samples)
            sol = DifferentialEquations.solve(jump_prob,DifferentialEquations.FunctionMap())
            max_time = sol.t[end]
            for i in 1:n_samples
                      outcome = round.(Int64,(sol(i*(max_time/(n_samples+1)))))
                      distribution_data[i,1] = data_scaling_RNAseq(outcome, c6)
            end
      return distribution_data
end
################################################




function generate_single_simulation_samples_m9_no_deact_no_k(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )

            c1 = 10^parameter_input[1]  #Activation
            c2 = 0.  #Deactivation
            c3 = 1.0                    #Expression
            c4 = 10^parameter_input[2]  #Degradation
            c5 = 10^parameter_input[3]     #Feedback baseline
            k =  3.    #Feedback k
            c6 = 0.977   #Scale factor
            G = 2
            P =(
                      (
                                sqrt(4*c1*c3*c4*G*k*(c1+c2+c5)+(c1*c4*k-c1*c3*G+c2*c4*k)^2)
                                +
                                (c1*c3*G-c1*c4*k-c2*c4*k)
                      )
                      /
                      (2*c4*(c1+c2+c5))
            )

            G1 = round(Int64,(

            (c1*G*(k+P))
            /
            (c1*(k+P)+c2*(k+P)+c5*P)
            ))

            P = round(Int64,(P))
            if P < 0
                      P = 0
            else
                      P = P
            end

            if G1 < 0
                      G1 = 0
            else
                      G1 = G1
            end
            G0 = 2 - G1
            starting_state = [G0,G1,P]
            time = 250000.0
            time_range= (0.0,time)
            #Define model
            #Set rate parameters
            rate1 = function (u,p,t)
                                                return (c1)*u[1]
            end

            affect! = function (integrator)
                      integrator.u[1] -= 1
                      integrator.u[2] += 1
            end
            jump1 = DifferentialEquations.ConstantRateJump(rate1,affect!)

            rate2 = function (u,p,t)
                      return (c2*u[2])+(u[2]*c5*((u[3]/k)/(1+(u[3]/k))))
            end

            affect! = function (integrator)
                      integrator.u[1] += 1
                      integrator.u[2] -= 1
            end
            jump2 = DifferentialEquations.ConstantRateJump(rate2,affect!)

            rate3 = function (u,p,t)
                      return (c3)*u[2]
             end
            affect! = function (integrator)
                      integrator.u[2] += 0
                      integrator.u[3] += 1
            end
            jump3 = DifferentialEquations.ConstantRateJump(rate3,affect!)

            rate4 = function (u,p,t)
                      return (c4)*u[3]
            end
            affect! = function (integrator)
                      integrator.u[3] -= 1
            end
            jump4 = DifferentialEquations.ConstantRateJump(rate4,affect!)

            prob = DifferentialEquations.DiscreteProblem(starting_state,time_range)
            jump_prob = DifferentialEquations.JumpProblem(prob,DifferentialEquations.Direct(),jump1,jump2,jump3,jump4)

            #Simulate multiple cells and take end point
            distribution_data = Array{Float64}(n_samples)
            sol = DifferentialEquations.solve(jump_prob,DifferentialEquations.FunctionMap())
            max_time = sol.t[end]
            for i in 1:n_samples
                      outcome = round.(Int64,(sol(i*(max_time/(n_samples+1)))))
                      distribution_data[i,1] = data_scaling_RNAseq(outcome, c6)
            end
      return distribution_data
end





#################################################

function generate_single_simulation_samples_m9_no_deact_no_k_save(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )

            c1 = 10^parameter_input[1]  #Activation
            c2 = 0.  #Deactivation
            c3 = 1.0                    #Expression
            c4 = 10^parameter_input[2]  #Degradation
            c5 = 10^parameter_input[3]     #Feedback baseline
            k =  3.    #Feedback k
            c6 = 0.977   #Scale factor
            G = 2
            P =(
                      (
                                sqrt(4*c1*c3*c4*G*k*(c1+c2+c5)+(c1*c4*k-c1*c3*G+c2*c4*k)^2)
                                +
                                (c1*c3*G-c1*c4*k-c2*c4*k)
                      )
                      /
                      (2*c4*(c1+c2+c5))
            )

            G1 = round(Int64,(

            (c1*G*(k+P))
            /
            (c1*(k+P)+c2*(k+P)+c5*P)
            ))

            P = round(Int64,(P))
            if P < 0
                      P = 0
            else
                      P = P
            end

            if G1 < 0
                      G1 = 0
            else
                      G1 = G1
            end
            G0 = 2 - G1
            starting_state = [G0,G1,P]
            time = 125000.0
            time_range= (0.0,time)
            #Define model
            #Set rate parameters
            rate1 = function (u,p,t)
                                                return (c1)*u[1]
            end

            affect! = function (integrator)
                      integrator.u[1] -= 1
                      integrator.u[2] += 1
            end
            jump1 = DifferentialEquations.ConstantRateJump(rate1,affect!)

            rate2 = function (u,p,t)
                      return (c2*u[2])+(u[2]*c5*((u[3]/k)/(1+(u[3]/k))))
            end

            affect! = function (integrator)
                      integrator.u[1] += 1
                      integrator.u[2] -= 1
            end
            jump2 = DifferentialEquations.ConstantRateJump(rate2,affect!)

            rate3 = function (u,p,t)
                      return (c3)*u[2]
             end
            affect! = function (integrator)
                      integrator.u[2] += 0
                      integrator.u[3] += 1
            end
            jump3 = DifferentialEquations.ConstantRateJump(rate3,affect!)

            rate4 = function (u,p,t)
                      return (c4)*u[3]
            end
            affect! = function (integrator)
                      integrator.u[3] -= 1
            end
            jump4 = DifferentialEquations.ConstantRateJump(rate4,affect!)

            prob = DifferentialEquations.DiscreteProblem(starting_state,time_range)
            jump_prob = DifferentialEquations.JumpProblem(prob,DifferentialEquations.Direct(),jump1,jump2,jump3,jump4)

            #Simulate multiple cells and take end point
            distribution_data = Array{Float64}(n_samples)
            sol = DifferentialEquations.solve(jump_prob,DifferentialEquations.FunctionMap())
            max_time = sol.t[end]
            for i in 1:n_samples
                      outcome = round.(Int64,(sol(i*(max_time/(n_samples+1)))))
                      distribution_data[i,1] = data_scaling_RNAseq(outcome, c6)
            end
      return distribution_data
end



################################################ m1_cell_cycle ##############################################################
function cell_cycle_simulation_generate_single_simulation(
    parameter_input::Array{Float64,1},
    cycle_time_vec::Array{Float64,1},
    cell_cycle_model
    )
    #define constants
    max_time = 40000.0
    time_1,time_2 = cycle_time_vec

    c1 = 10^parameter_input[1]
    c2 = 10^parameter_input[2]
    c3 = 1.0
    c4 = 10^parameter_input[3]
    c5 = 1.
    p = [c1,c2,c3,c4]

    #storage of output
    solution = []
    time = []

    #Define initial starting state (starting in 2-state system)
    G = 2
    G1 = round(Int64,(2*c1)/(c1+c2))
    G0 = G - G1
    P = round(Int64,((c3/c4)*((2*c1)/(c1+c2))))
    starting_state = [G0,G1,P]

    time_range = (0.0,time_1)
    end_time = 0
    while end_time <= max_time
        #Cell cycle part 1
        prob = DiscreteProblem(starting_state,time_range,p)
        jump_prob = JumpProblem(prob,Direct(),cell_cycle_model)
        sol = solve(jump_prob,FunctionMap())
        end_time = sol.t[end]
        solution = vcat(solution,sol.u)
        time = vcat(time,sol.t)

        #Cell cycle part 2
        starting_state = sol.u[end] #Redefine the starting state based on end point of previous simulation
        starting_state[1] = starting_state[1]*2 #Double number of genes, keeping their current status of activitity (on/off)
        starting_state[2] = starting_state[2]*2
        time_range= (end_time,end_time+time_2)

        prob = DiscreteProblem(starting_state,time_range,p)
        jump_prob = JumpProblem(prob,Direct(),cell_cycle_model)
        sol = solve(jump_prob,FunctionMap())
        end_time = sol.t[end]
        solution = vcat(solution,sol.u)
        time = vcat(time,sol.t)

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
        time_range=(end_time,end_time+time_1)
    end
    return(solution,time)
end
########
#Simulated product throughout time and randomly samples to produce a distribution
function cell_cycle_simulation_generate_single_simulation_samples_m1_only_one(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )
            max_time = 250000.0
            cycle_time_vec = [83333.,166666.]

            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            time_1,time_2 = cycle_time_vec
            c1 = 10^parameter_input[1]
            c2 = 10^parameter_input[2]
            c3 = 1.0
            c4 = 10^parameter_input[3]
            c5 = 1.
            p = [c1,c2,c3,c4]
            #Define initial starting state (starting in 2-state system)
            G = 2
            G1 = round(Int64,(2*c1)/(c1+c2))
            G0 = G - G1
            P = round(Int64,((c3/c4)*((2*c1)/(c1+c2))))
            starting_state = [G0,G1,P]
            time_range = (0.0,time_1)
            end_time = 0
            output_samples = []
            start_time = 10
            end_time = 2000
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            while end_time <= max_time
                       #Cell cycle part 1
                       prob = DiscreteProblem(starting_state,time_range,p)
                       jump_prob = JumpProblem(prob,Direct(),base_model)
                       sol = solve(jump_prob,FunctionMap())
                       start_time = sol.t[1]
                       end_time = sol.t[end]
                       sub_sample_points = sample_points[[x >= start_time for x in sample_points]]
                       sub_sample_points = sub_sample_points[[x < end_time for x in sub_sample_points]]
                       if length(sub_sample_points) >= 1
                           output = sol(sub_sample_points)
                           # print(output)
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
                       time_range=(end_time,end_time+time_1)
           end
      return (float([i for i in output_samples]))
end



########
#Simulated product throughout time and randomly samples to produce a distribution
function cell_cycle_simulation_generate_single_simulation_samples_m1_8333_16666(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )
            cycle_time_vec= [8333.,16666.]
            max_time = 250000.0
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            time_1,time_2 = cycle_time_vec
            c1 = 10^parameter_input[1]
            c2 = 10^parameter_input[2]
            c3 = 1.0
            c4 = 10^parameter_input[3]
            c5 = 1.
            p = [c1,c2,c3,c4]
            #Define initial starting state (starting in 2-state system)
            G = 2
            G1 = round(Int64,(2*c1)/(c1+c2))
            G0 = G - G1
            P = round(Int64,((c3/c4)*((2*c1)/(c1+c2))))
            starting_state = [G0,G1,P]
            time_range = (0.0,time_1)
            end_time = 0
            output_samples = []
            start_time = 10
            end_time = 2000
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            while end_time <= max_time
                       #Cell cycle part 1
                       prob = DiscreteProblem(starting_state,time_range,p)
                       jump_prob = JumpProblem(prob,Direct(),base_model)
                       sol = solve(jump_prob,FunctionMap())
                       start_time = sol.t[1]
                       end_time = sol.t[end]
                       sub_sample_points = sample_points[[x >= start_time for x in sample_points]]
                       sub_sample_points = sub_sample_points[[x < end_time for x in sub_sample_points]]
                       if length(sub_sample_points) >= 1
                           output = sol(sub_sample_points)
                           # print(output)
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
                       time_range=(end_time,end_time+time_1)
           end
      return (float([i for i in output_samples]))
end
######################################################################################################################
#Simulated product throughout time and randomly samples to produce a distribution
function cell_cycle_simulation_generate_single_simulation_samples_m1_16666_33332(
            parameter_input::Vector{Float64},
            n_samples::Int64,
            data_scaling_function::Function
            )
            cycle_time_vec= [16666.,33332.]
            max_time = 250000.0
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            time_1,time_2 = cycle_time_vec
            c1 = 10^parameter_input[1]
            c2 = 10^parameter_input[2]
            c3 = 1.0
            c4 = 10^parameter_input[3]
            c5 = 1.
            p = [c1,c2,c3,c4]
            #Define initial starting state (starting in 2-state system)
            G = 2
            G1 = round(Int64,(2*c1)/(c1+c2))
            G0 = G - G1
            P = round(Int64,((c3/c4)*((2*c1)/(c1+c2))))
            starting_state = [G0,G1,P]
            time_range = (0.0,time_1)
            end_time = 0
            output_samples = []
            start_time = 10
            end_time = 2000
            sample_points = StatsBase.sample(0.0:max_time,n_samples)
            while end_time <= max_time
                       #Cell cycle part 1
                       prob = DiscreteProblem(starting_state,time_range,p)
                       jump_prob = JumpProblem(prob,Direct(),base_model)
                       sol = solve(jump_prob,FunctionMap())
                       start_time = sol.t[1]
                       end_time = sol.t[end]
                       sub_sample_points = sample_points[[x >= start_time for x in sample_points]]
                       sub_sample_points = sub_sample_points[[x < end_time for x in sub_sample_points]]
                       if length(sub_sample_points) >= 1
                           output = sol(sub_sample_points)
                           # print(output)
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
                       time_range=(end_time,end_time+time_1)
           end
      return (float([i for i in output_samples]))
end
