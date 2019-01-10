using PyCall, PyPlot; @pyimport seaborn as sns
function plot_one_abacus_file_m1(pathtofile_lisi, pathtoplot, number_abc_smc, number_particles, number_unknown_params)
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
                        end
                    end
                end
            end
        end
    end
    a=a[1:end-2]
    b=b[1:end-2]
    c=c[1:end-2]
    # d=d[1:end-2]
    a2=a2[1:end-2]
    b2=b2[1:end-2]
    c2=c2[1:end-2]
    # d2=d2[1:end-2]
    a3=a3[1:end-2]
    b3=b3[1:end-2]
    c3=c3[1:end-2]
    # d3=d3[1:end-2]
    PyPlot.ioff()
    fig, ax = PyPlot.subplots(figsize=(15,15),ncols=3, nrows=3)
    left   =  0.125  # the left side of the subplots of the figure
    right  =  0.9    # the right side of the subplots of the figure
    bottom =  0.1    # the bottom of the subplots of the figure
    top    =  0.9    # the top of the subplots of the figure
    wspace =  .5     # the amount of width reserved for blank space between subplots
    hspace =  1.1    # the amount of height reserved for white space between subplots
    PyPlot.subplots_adjust(
        left    =  left,
        bottom  =  bottom,
        right   =  right,
        top     =  top,
        wspace  =  wspace,
        hspace  =  hspace
    )
    y_title_margin = 1.2
    PyPlot.suptitle("Original vs Normalized vs Standardized", y = 1.09, fontsize=20)
    sns.distplot(a,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c act",ax=ax[1,1])
    sns.distplot(b,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c de-act",ax=ax[2,1])
    sns.distplot(c,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c deg",ax=ax[3,1])
    # sns.distplot(d,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "scaling factor",ax=ax[4,1])

    sns.distplot(a2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c act",ax=ax[1,2])
    sns.distplot(b2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c de-act",ax=ax[2,2])
    sns.distplot(c2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c deg",ax=ax[3,2])
    # sns.distplot(d2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "scaling factor",ax=ax[4,2])

    sns.distplot(a3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "c act",ax=ax[1,3])
    sns.distplot(b3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "c de-act",ax=ax[2,3])
    sns.distplot(c3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "c deg",ax=ax[3,3])
    # sns.distplot(d3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "scaling factor",ax=ax[4,3])

    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[1,1])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[1,2])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[1,3])


    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[2,1])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[2,2])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[2,3])


    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[3,1])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[3,2])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[3,3])
    #
    #
    # sns.distplot([-2,22],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[4,1])
    # sns.distplot([-2,22],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[4,2])
    # sns.distplot([-2,22],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[4,3])



    PyPlot.savefig(pathtoplot)
    return  a,b,c,a2,b2,c2,a3,b3,c3
end




function plot_one_abacus_file_m8(pathtofile_lisi, pathtoplot, number_abc_smc, number_particles, number_unknown_params)
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
                        if count_particles[1]<=100
                            ln_arr =split(ln)
                            push!(a,parse(Float64,ln_arr[1]))
                            push!(b,parse(Float64,ln_arr[2]))
                            push!(c,parse(Float64,ln_arr[3]))
                            push!(d,parse(Float64,ln_arr[4]))
                            push!(e,parse(Float64,ln_arr[5]))
                        end
                    elseif (counter_rejection_abc[1]==2)
                        count_particles[1]+=1
                        if count_particles[1]<=100
                            ln_arr =split(ln)
                            push!(a2,parse(Float64,ln_arr[1]))
                            push!(b2,parse(Float64,ln_arr[2]))
                            push!(c2,parse(Float64,ln_arr[3]))
                            push!(d2,parse(Float64,ln_arr[4]))
                            push!(e2,parse(Float64,ln_arr[5]))
                        end
                    elseif (counter_rejection_abc[1]==3)
                        count_particles[1]+=1
                        if count_particles[1]<=100
                            ln_arr =split(ln)
                            push!(a3,parse(Float64,ln_arr[1]))
                            push!(b3,parse(Float64,ln_arr[2]))
                            push!(c3,parse(Float64,ln_arr[3]))
                            push!(d3,parse(Float64,ln_arr[4]))
                            push!(e3,parse(Float64,ln_arr[5]))
                        end
                    end
                end
            end
        end
    end
    a=a[1:end-2]
    b=b[1:end-2]
    c=c[1:end-2]
    d=d[1:end-2]
    e=e[1:end-2]
    a2=a2[1:end-2]
    b2=b2[1:end-2]
    c2=c2[1:end-2]
    d2=d2[1:end-2]
    e2=e2[1:end-2]
    a3=a3[1:end-2]
    b3=b3[1:end-2]
    c3=c3[1:end-2]
    d3=d3[1:end-2]
    e3=e3[1:end-2]
    PyPlot.ioff()
    fig, ax = PyPlot.subplots(figsize=(15,15),ncols=3, nrows=5)
    left   =  0.125  # the left side of the subplots of the figure
    right  =  0.9    # the right side of the subplots of the figure
    bottom =  0.1    # the bottom of the subplots of the figure
    top    =  0.9    # the top of the subplots of the figure
    wspace =  .5     # the amount of width reserved for blank space between subplots
    hspace =  1.1    # the amount of height reserved for white space between subplots
    PyPlot.subplots_adjust(
        left    =  left,
        bottom  =  bottom,
        right   =  right,
        top     =  top,
        wspace  =  wspace,
        hspace  =  hspace
    )
    y_title_margin = 1.2
    PyPlot.suptitle("Original vs Normalized vs Standardized", y = 1.09, fontsize=20)
    sns.distplot(a,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c act",ax=ax[1,1])
    sns.distplot(b,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c de-act",ax=ax[2,1])
    sns.distplot(c,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c deg",ax=ax[3,1])
    sns.distplot(d,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "feedback",ax=ax[4,1])
    sns.distplot(e,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "scaling factor",ax=ax[5,1])

    sns.distplot(a2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c act",ax=ax[1,2])
    sns.distplot(b2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c de-act",ax=ax[2,2])
    sns.distplot(c2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c deg",ax=ax[3,2])
    sns.distplot(d2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "feedback",ax=ax[4,2])
    sns.distplot(e2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "scaling factor",ax=ax[5,2])

    sns.distplot(a3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "c act",ax=ax[1,3])
    sns.distplot(b3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "c de-act",ax=ax[2,3])
    sns.distplot(c3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "c deg",ax=ax[3,3])
    sns.distplot(d3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "feedback",ax=ax[4,3])
    sns.distplot(e3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "scaling factor",ax=ax[5,3])

    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[1,1])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[1,2])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[1,3])


    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[2,1])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[2,2])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[2,3])


    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[3,1])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[3,2])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[3,3])

    sns.distplot([-4,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[4,1])
    sns.distplot([-4,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[4,2])
    sns.distplot([-4,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[4,3])

    sns.distplot([-2,12],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[5,1])
    sns.distplot([-2,12],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[5,2])
    sns.distplot([-2,12],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[5,3])
    PyPlot.savefig(pathtoplot)

end




function plot_one_abacus_file_m9(pathtofile_lisi, pathtoplot, number_abc_smc, number_particles, number_unknown_params)
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
                        if count_particles[1]<=100
                            ln_arr =split(ln)
                            push!(a,parse(Float64,ln_arr[1]))
                            push!(b,parse(Float64,ln_arr[2]))
                            push!(c,parse(Float64,ln_arr[3]))
                            push!(d,parse(Float64,ln_arr[4]))
                            push!(e,parse(Float64,ln_arr[5]))
                            push!(f,parse(Float64,ln_arr[6]))
                        end
                    elseif (counter_rejection_abc[1]==2)
                        count_particles[1]+=1
                        if count_particles[1]<=100
                            ln_arr =split(ln)
                            push!(a2,parse(Float64,ln_arr[1]))
                            push!(b2,parse(Float64,ln_arr[2]))
                            push!(c2,parse(Float64,ln_arr[3]))
                            push!(d2,parse(Float64,ln_arr[4]))
                            push!(e2,parse(Float64,ln_arr[5]))
                            push!(f2,parse(Float64,ln_arr[6]))

                        end
                    elseif (counter_rejection_abc[1]==3)
                        count_particles[1]+=1
                        if count_particles[1]<=100
                            ln_arr =split(ln)
                            push!(a3,parse(Float64,ln_arr[1]))
                            push!(b3,parse(Float64,ln_arr[2]))
                            push!(c3,parse(Float64,ln_arr[3]))
                            push!(d3,parse(Float64,ln_arr[4]))
                            push!(e3,parse(Float64,ln_arr[5]))
                            push!(f3,parse(Float64,ln_arr[6]))

                        end
                    end
                end
            end
        end
    end
    a=a[1:end-2]
    b=b[1:end-2]
    c=c[1:end-2]
    d=d[1:end-2]
    e=e[1:end-2]
    f=f[1:end-2]

    a2=a2[1:end-2]
    b2=b2[1:end-2]
    c2=c2[1:end-2]
    d2=d2[1:end-2]
    e2=e2[1:end-2]
    f2=f2[1:end-2]
    a3=a3[1:end-2]
    b3=b3[1:end-2]
    c3=c3[1:end-2]
    d3=d3[1:end-2]
    e3=e3[1:end-2]
    f3=f3[1:end-2]

    PyPlot.ioff()
    fig, ax = PyPlot.subplots(figsize=(15,15),ncols=3, nrows=6)
    left   =  0.125  # the left side of the subplots of the figure
    right  =  0.9    # the right side of the subplots of the figure
    bottom =  0.1    # the bottom of the subplots of the figure
    top    =  0.9    # the top of the subplots of the figure
    wspace =  .5     # the amount of width reserved for blank space between subplots
    hspace =  1.1    # the amount of height reserved for white space between subplots
    PyPlot.subplots_adjust(
        left    =  left,
        bottom  =  bottom,
        right   =  right,
        top     =  top,
        wspace  =  wspace,
        hspace  =  hspace
    )
    y_title_margin = 1.2
    PyPlot.suptitle("Original vs Normalized vs Standardized", y = 1.09, fontsize=20)
    sns.distplot(a,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c act",ax=ax[1,1])
    sns.distplot(b,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c de-act",ax=ax[2,1])
    sns.distplot(c,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "c deg",ax=ax[3,1])
    sns.distplot(d,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "feedback baseline",ax=ax[4,1])
    sns.distplot(e,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "feedback k",ax=ax[5,1])
    sns.distplot(f,bins=20, kde= true, hist=true,rug= true ,  color = "thistle", label = "lisi",axlabel = "scaling factor",ax=ax[6,1])

    sns.distplot(a2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c act",ax=ax[1,2])
    sns.distplot(b2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c de-act",ax=ax[2,2])
    sns.distplot(c2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "c deg",ax=ax[3,2])
    sns.distplot(d2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "feedback baseline",ax=ax[4,2])
    sns.distplot(e2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "feedback k",ax=ax[5,2])
    sns.distplot(f2,bins=20, kde= true, hist=true,rug= true ,  color = "orchid", label = "lisi",axlabel = "scaling factor",ax=ax[6,2])


    sns.distplot(a3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "c act",ax=ax[1,3])
    sns.distplot(b3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "c de-act",ax=ax[2,3])
    sns.distplot(c3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "c deg",ax=ax[3,3])
    sns.distplot(d3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "feedback baseline",ax=ax[4,3])
    sns.distplot(e3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "feedback k",ax=ax[5,3])
    sns.distplot(f3,bins=20, kde= true, hist=true,rug= true ,  color = "crimson", label = "lisi",axlabel = "scaling factor",ax=ax[6,3])


    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[1,1])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[1,2])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[1,3])



    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[2,1])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[2,2])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[2,3])


    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[3,1])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[3,2])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[3,3])


    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[4,1])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[4,2])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[4,3])


    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[5,1])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[5,2])
    sns.distplot([-7,1],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[5,3])

    sns.distplot([-2,12],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[5,1])
    sns.distplot([-2,12],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[5,2])
    sns.distplot([-2,12],bins=20,  kde= false, hist=true,rug= false ,  color = "white", label = "lisi",ax=ax[5,3])
    PyPlot.savefig(pathtoplot)

end
