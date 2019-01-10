mutable struct hazard_fuction
    k_i::Float64
    current_state::Array{Float64,1}
    new_state::Array{Float64,1}
end

function gill_li(
    pre_reaction::Array{Float64,2}, post_reaction::Array{Float64,2},
    k::Array{Float64,1}, u0::Array{Float64,1}, tmax::Float64)
    if size(pre_reaction)!=size(post_reaction)
        print("Dim mismatch.")
        return
    end

    number_species = size(pre_reaction)[1]
    number_reactions = size(pre_reaction)[2]
    #init temp with start point
    temp_state = u0
    temp_time = 0.

    S =  post_reaction - pre_reaction
    hazard_functions = Dict{Float64,hazard_fuction}
    for reaction_i in 1:number_reactions
            k_i = k[reaction_i]
            for
            hazard_fuction = [k_i, temp_state, ]
end


pre_reaction = [1 2 3;1 2 3;1 2 3;1. 2 3]
post_reaction = [1 2 0;1 0 3;1 0 3;1. 0 3]

ktest = [1.,2,2,1]
u0test = [1.,1.,1.]
tmaxtest =200.
gill_li(pre_reaction,post_reaction,ktest,u0test,tmaxtest)
number_species=2


x=[1,2,2]
find(isequal(2),x)
inds=[]
for i in 1:length(x)
    if(x[i]!=2)
        push!(inds,i)
    end
end

hazard_functions = Dict{Float64,hazard_fuction}

for reaction in
