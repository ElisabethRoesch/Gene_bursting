function read_one_m8(pathtofile_lisi,  number_abc_smc, number_particles, number_unknown_params)
    counter_headerlines = [0]
    counter_rejection_abc=[0]
    count_particles = [0]
    #enter number of accepted particles and number of unknown parameters
    final_abc =Array{Float64,2}(number_particles,number_unknown_params)
    a=[]
    b=[]
    c=[]
    d=[]
    e=[]

    a2=[]
    b2=[]
    c2=[]
    d2=[]
    e2=[]


    a3=[]
    b3=[]
    c3=[]
    d3=[]
    e3=[]

    di1=[]
    di2=[]
    di3=[]

    w1=[]
    w2=[]
    w3=[]
    open(pathtofile_lisi) do file
        for ln in eachline(file)
            if startswith(ln, "-----")
                counter_headerlines[1]=1
                count_particles[1]=0
                counter_rejection_abc[1]+=1
                # println("$(length(ln)), $(ln)")
            elseif counter_headerlines[1]>0
                counter_headerlines[1]+=1
                if counter_headerlines[1] >= 5
                    if (counter_rejection_abc[1]==1)
                        count_particles[1]+=1
                        if count_particles[1]<=number_particles
                            ln_arr =split(ln)
                            push!(a,parse(Float64,ln_arr[1]))
                            push!(b,parse(Float64,ln_arr[2]))
                            push!(c,parse(Float64,ln_arr[3]))
                            push!(d,parse(Float64,ln_arr[4]))
                            # push!(e,parse(Float64,ln_arr[5]))
                        end
                    elseif (counter_rejection_abc[1]==2)
                        count_particles[1]+=1
                        if count_particles[1]<=number_particles
                            ln_arr =split(ln)
                            push!(a2,parse(Float64,ln_arr[1]))
                            push!(b2,parse(Float64,ln_arr[2]))
                            push!(c2,parse(Float64,ln_arr[3]))
                            push!(d2,parse(Float64,ln_arr[4]))
                            # push!(e2,parse(Float64,ln_arr[5]))
                        end
                    elseif (counter_rejection_abc[1]==3)
                        count_particles[1]+=1
                        if count_particles[1]<=number_particles
                            ln_arr =split(ln)
                            push!(a3,parse(Float64,ln_arr[1]))
                            push!(b3,parse(Float64,ln_arr[2]))
                            push!(c3,parse(Float64,ln_arr[3]))
                            push!(d3,parse(Float64,ln_arr[4]))
                            # push!(e3,parse(Float64,ln_arr[5]))
                        elseif count_particles[1]==number_particles+1
                            ln_arr =split(ln)
                            for i in ln_arr
                                push!(di3,parse(Float64,i))
                            end
                            # println("di3",di3)
                        elseif count_particles[1]==number_particles+2
                            ln_arr =split(ln)
                            for i in ln_arr
                                push!(w3,parse(Float64,i))
                            end
                            # println("w3",w3)
                        end
                    end
                end
            end
        end
    end
    return a3,b3,c3,d3,e3,di3,w3
end


function read_one_m9(pathtofile_lisi,  number_abc_smc, number_particles, number_unknown_params)
    counter_headerlines = [0]
    counter_rejection_abc=[0]
    count_particles = [0]
    #enter number of accepted particles and number of unknown parameters
    final_abc =Array{Float64,2}(number_particles,number_unknown_params)
    a=[]
    b=[]
    c=[]
    d=[]
    e=[]
    f=[]


    a2=[]
    b2=[]
    c2=[]
    d2=[]
    e2=[]
    f2=[]


    a3=[]
    b3=[]
    c3=[]
    d3=[]
    e3=[]
    f3=[]

    di1=[]
    di2=[]
    di3=[]

    w1=[]
    w2=[]
    w3=[]
    open(pathtofile_lisi) do file
        for ln in eachline(file)
            if startswith(ln, "-----")
                counter_headerlines[1]=1
                count_particles[1]=0
                counter_rejection_abc[1]+=1
                # println("$(length(ln)), $(ln)")
            elseif counter_headerlines[1]>0
                counter_headerlines[1]+=1
                if counter_headerlines[1] >= 5
                    if (counter_rejection_abc[1]==1)
                        count_particles[1]+=1
                        if count_particles[1]<=number_particles
                            ln_arr =split(ln)
                            push!(a,parse(Float64,ln_arr[1]))
                            push!(b,parse(Float64,ln_arr[2]))
                            push!(c,parse(Float64,ln_arr[3]))
                            push!(d,parse(Float64,ln_arr[4]))
                            push!(e,parse(Float64,ln_arr[5]))
                            # push!(f,parse(Float64,ln_arr[6]))

                        end
                    elseif (counter_rejection_abc[1]==2)
                        count_particles[1]+=1
                        if count_particles[1]<=number_particles
                            ln_arr =split(ln)
                            push!(a2,parse(Float64,ln_arr[1]))
                            push!(b2,parse(Float64,ln_arr[2]))
                            push!(c2,parse(Float64,ln_arr[3]))
                            push!(d2,parse(Float64,ln_arr[4]))
                            push!(e2,parse(Float64,ln_arr[5]))
                            # push!(f2,parse(Float64,ln_arr[6]))

                        end
                    elseif (counter_rejection_abc[1]==3)
                        count_particles[1]+=1
                        if count_particles[1]<=number_particles
                            ln_arr =split(ln)
                            push!(a3,parse(Float64,ln_arr[1]))
                            push!(b3,parse(Float64,ln_arr[2]))
                            push!(c3,parse(Float64,ln_arr[3]))
                            push!(d3,parse(Float64,ln_arr[4]))
                            push!(e3,parse(Float64,ln_arr[5]))
                            # push!(f3,parse(Float64,ln_arr[6]))

                        elseif count_particles[1]==number_particles+1
                            ln_arr =split(ln)
                            for i in ln_arr
                                push!(di3,parse(Float64,i))
                            end
                            # println("di3",di3)
                        elseif count_particles[1]==number_particles+2
                            ln_arr =split(ln)
                            for i in ln_arr
                                push!(w3,parse(Float64,i))
                            end
                            # println("w3",w3)
                        end
                    end
                end
            end
        end
    end
    return a3,b3,c3,d3,e3,f3,di3,w3
end





function read_one_m1(pathtofile_lisi,  number_abc_smc, number_particles, number_unknown_params)
    counter_headerlines = [0]
    counter_rejection_abc=[0]
    count_particles = [0]
    #enter number of accepted particles and number of unknown parameters
    final_abc =Array{Float64,2}(number_particles,number_unknown_params)
    a=[]
    b=[]
    c=[]
    d=[]

    a2=[]
    b2=[]
    c2=[]
    d2=[]


    a3=[]
    b3=[]
    c3=[]
    d3=[]

    di1=[]
    di2=[]
    di3=[]

    w1=[]
    w2=[]
    w3=[]
    open(pathtofile_lisi) do file
        for ln in eachline(file)
            if startswith(ln, "-----")
                counter_headerlines[1]=1
                count_particles[1]=0
                counter_rejection_abc[1]+=1
                # println("$(length(ln)), $(ln)")
            elseif counter_headerlines[1]>0
                counter_headerlines[1]+=1
                if counter_headerlines[1] >= 5
                    if (counter_rejection_abc[1]==1)
                        count_particles[1]+=1
                        if count_particles[1]<=number_particles
                            ln_arr =split(ln)
                            push!(a,parse(Float64,ln_arr[1]))
                            push!(b,parse(Float64,ln_arr[2]))
                            push!(c,parse(Float64,ln_arr[3]))
                            # push!(d,parse(Float64,ln_arr[4]))
                        end
                    elseif (counter_rejection_abc[1]==2)
                        count_particles[1]+=1
                        if count_particles[1]<=number_particles
                            ln_arr =split(ln)
                            push!(a2,parse(Float64,ln_arr[1]))
                            push!(b2,parse(Float64,ln_arr[2]))
                            push!(c2,parse(Float64,ln_arr[3]))
                            # push!(d2,parse(Float64,ln_arr[4]))
                        end
                    elseif (counter_rejection_abc[1]==3)
                        count_particles[1]+=1
                        if count_particles[1]<=number_particles
                            ln_arr =split(ln)
                            push!(a3,parse(Float64,ln_arr[1]))
                            push!(b3,parse(Float64,ln_arr[2]))
                            push!(c3,parse(Float64,ln_arr[3]))
                            # push!(d3,parse(Float64,ln_arr[4]))
                        elseif count_particles[1]==number_particles+1
                            ln_arr =split(ln)
                            for i in ln_arr
                                push!(di3,parse(Float64,i))
                            end
                            # println("di3",di3)
                        elseif count_particles[1]==number_particles+2
                            ln_arr =split(ln)
                            for i in ln_arr
                                push!(w3,parse(Float64,i))
                            end
                            # println("w3",w3)
                        end
                    end
                end
            end
        end
    end
    # print("dist")
    # print(di3)
    return a3,b3,c3,d3,di3,w3
end
