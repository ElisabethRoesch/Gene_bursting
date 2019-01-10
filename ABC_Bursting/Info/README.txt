________________________________________________________________________________
ABaCus package
---------------

  ABaCus.jl
  io.jl

Info: - A copy of the package compatible with all other code used is included in
        the Code folder.
        Any updates made to the ABaCus package post June 2017 may be incompatible

________________________________________________________________________________
Gillespie simulations of expressed product distributions:
----------------------------------------------------------

  Gillespie_model_1_fixedprod.jl  standard model
  Gillespie_model_1_fixedprod_cellcycle.jl  standard model with cell cycle
  Gillespie_model_8_fixedprod.jl  negative feedback regulated by activation rate
  Gillespie_model_9_fixedprod.jl  positive feedback regulated by de-activation rate
  data_scaling.jl includes functions to scale absolute product number

Info: - See transcriptional_bursting_model.pdf for more details of models and diagrams
      - Gillespie_model_8_fixedprod.jl and Gillespie_model_9_fixedprod.jl will
        need modifying for use with the updated DifferentialEquations.jl package


________________________________________________________________________________
Running ABC-SMC
-----------------

  ABC_SMC_RNAseq_fixedparams.jl
    Acts as a template form to fill in details of the run, e.g. specify genes, data, priors

  run_singleABC.sh
    Is a template for performing a single run of ABC on one CPU

  run_ABC.sh
    Is an interface to allow running the ABC on several servers, specified by args
    Takes in arguments to allow:
      Defining the number of processes to run
      Running of the same ABC_SMC_RNAseq_fixedparams.jl across different servers
      Indexing the full gene name list to run across different servers
      Selecting output file name
      Selecting log file locations
    This file should not be edited

  run_singleABC.sh
    Is a template for performing a single run of ABC on several servers using
      run_ABC.sh
