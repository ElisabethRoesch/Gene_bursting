    priors = [Distributions.Uniform(-3.0,3.0),
        Distributions.Uniform(-3.0,3.0),
        Distributions.Uniform(-3.0,3.0),
        Distributions.Uniform(-3.0,0.5),
        Distributions.Uniform(0.0,5.0)] 

    #Define parameters for SMC
    n_params = 5 
    n_particles = 100
    threshold = [5.,1.5,0.8]
