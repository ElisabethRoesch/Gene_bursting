using PyPlot
PyPlot.ion()

serum_act =[(0.304, 8.1206e-21), (0.192, 1.46909e-8), (0.232, 2.65981e-12), (0.128, 0.000485207), (0.096, 0.0185154), (0.08, 0.0774203), (0.146, 3.95994e-5), (0.156, 8.5375e-6), (0.2, 2.98603e-9), (0.172, 5.93428e-7), (0.108, 0.00533811)]
serum_deact = [(0.782, 2.34014e-135), (0.6, 7.37239e-80), (0.626, 6.82684e-87), (0.456, 2.632e-46), (0.42, 2.39141e-39), (0.3, 2.77124e-20), (0.366, 5.54439e-30), (0.49, 2.10548e-53), (0.38, 2.75044e-32), (0.3, 2.77124e-20), (0.152, 1.5965e-5)]
serum_deg =[(0.512, 2.8799e-58), (0.41, 1.62194e-37), (0.45, 4.16611e-45), (0.228, 6.77398e-12), (0.292, 3.07373e-19), (0.342, 3.11251e-26), (0.154, 1.16986e-5), (0.374, 2.73915e-31), (0.132, 0.000286057), (0.182, 9.82364e-8), (0.38, 2.75044e-32)]

i_act = [(0.084, 0.0554786), (0.142, 7.1101e-5), (0.388, 1.21266e-33), (0.348, 3.79863e-27), (0.33, 1.87249e-24), (0.342, 3.11251e-26)]
i_deact = [(0.15, 2.16991e-5), (0.08, 0.0774203), (0.488, 5.68773e-53), (0.306, 4.36915e-21), (0.52, 4.34214e-60), (0.194, 9.92446e-9)]
i_deg = [(0.224, 1.69737e-11), (0.062, 0.282889), (0.136, 0.000165927), (0.292, 3.07373e-19), (0.176, 2.92567e-7), (0.216, 1.01498e-10)]

a= reinterpret(Float64, serum_act, (2, 11))
b= reinterpret(Float64, serum_deact, (2, 11))
c= reinterpret(Float64, serum_deg, (2, 11))
d= reinterpret(Float64, i_act, (2,6))
e= reinterpret(Float64, i_deact, (2,6))
f= reinterpret(Float64, i_deg,(2,6))

mean(a[1,:])
mean(b[1,:])
mean(c[1,:])

fig, ax = PyPlot.subplots(figsize=(3,10))
  left   =  0.3  # the left side of the subplots of the figure
  right  =  0.9    # the right side of the subplots of the figure
  bottom =  0.2    # the bottom of the subplots of the figure
  top    =  0.9    # the top of the subplots of the figure
  wspace =  .5     # the amount of width reserved for blank space between subplots
  hspace =  .5    # the amount of height reserved for white space between subplots
  PyPlot.subplots_adjust(
      left    =  left,
      bottom  =  bottom,
      right   =  right,
      top     =  top,
      wspace  =  wspace,
      hspace  =  hspace
  )
    # ax[:set_xlim]([0,0.8])
    ax[:tick_params](labelsize=19)
    ax[:set_ylabel]("KS distance",fontsize=19)
    PyPlot.scatter(fill(0.,length(a[1,:])),a[1,:],alpha=0.5, s =300,color ="grey")
    PyPlot.scatter(fill(1.,length(a[1,:])),b[1,:],alpha=0.5, s =300,color ="grey")
    PyPlot.scatter(fill(2.,length(a[1,:])),c[1,:],alpha=0.5,s =300,color ="grey")
    PyPlot.axhline([0.00544],color ="orange",linewidth=3, linestyle="--",label = "Critical value")
    ax[:set_xlim]([-0.4,2.4])
    ax[:set_ylim]([-0.05,0.85])

    ax[:legend](fontsize = 18,loc="top right",bbox_to_anchor=(1.05, 1.15))
    My = matplotlib[:ticker][:MultipleLocator](1) # Define interval of major ticks
    ax[:xaxis][:set_major_locator](My)
    # ax[:spines]["right"][:set_color]("none")
    # ax[:spines]["top"][:set_color]("none")

    ax[:set_xticklabels](["no","Activation","Deactivation","Degradation"],rotation=90)
    savefig("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/params/serumgenes_ksh.pdf")


fig, ax = PyPlot.subplots(figsize=(3,10))
      left   =  0.3  # the left side of the subplots of the figure
      right  =  0.9    # the right side of the subplots of the figure
      bottom =  0.2    # the bottom of the subplots of the figure
      top    =  0.9    # the top of the subplots of the figure
      wspace =  .5     # the amount of width reserved for blank space between subplots
      hspace =  .5    # the amount of height reserved for white space between subplots
      PyPlot.subplots_adjust(
          left    =  left,
          bottom  =  bottom,
          right   =  right,
          top     =  top,
          wspace  =  wspace,
          hspace  =  hspace
      )
    ax[:tick_params](labelsize=19)
    ax[:set_ylabel]("KS distance",fontsize=19)
    ax[:set_xlim]([-0.4,2.4])
    ax[:set_ylim]([-0.05,0.85])
    ax[:scatter](fill(0,length(d[1,:])),d[1,:],alpha=0.5,  s =300,color ="grey")
    ax[:scatter](fill(1,length(d[1,:])),e[1,:],alpha=0.5,  s =300,color ="grey")
    ax[:scatter](fill(2,length(d[1,:])),f[1,:],alpha=0.5,  s =300,color ="grey")
    ax[:axhline]([0.00544],color ="orange",linewidth=3, linestyle="--",label = "Critical value")
    # ax[:legend](fontsize = 18,loc="top right",bbox_to_anchor=(1.05, 1.15))
    My = matplotlib[:ticker][:MultipleLocator](1) # Define interval of major ticks
    ax[:xaxis][:set_major_locator](My)
    ax[:set_xticklabels](["no","Activation","Deactivation","Degradation"],rotation=90)

    savefig("/project/home17/er4517/project_3/gene_bursting/ABC_Bursting/plots/params/2igenes_ksh.pdf")

# Axes.set_yticks
# bbox_to_anchor=(1.05, 1), loc=2,
    # ax[2,1][:scatter]([truth[1]], [truth[2]],color ="gold", marker = "x",linewidth=3,  s =600, zorder =120)
    # My = matplotlib[:ticker][:MultipleLocator](1) # Define interval of major ticks
    # ax[2,1][:yaxis][:set_major_locator](My)
    # ax[1,1][:yaxis][:set_major_locator](My)
    # ax[1,1][:set_yticklabels](["Activation","Activation","Deactivation","Degradation"])
    # ax[2,1][:set_yticklabels](["Activation","Activation","Deactivation","Degradation"])
