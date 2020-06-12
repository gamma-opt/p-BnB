"""
    MIP_initial_parameters
Stores attributes for generating MIP JuMP model. Has the following fields:
* `num_scen::Int`:                              Number of scenarios used in the model
* `scen_prob::Vector{Float64}`:                 Vector of probabilities associated with each scenario
* `num_int_var::Int`:                           Number of integer variables associated with each scenario
* `num_cont_var::Int`:                          Number of continuous variables associated with each scenario
* `num_cont_var::Int`:                          Number of constraints associated with each scenario
* `quad_mat_dens::Float64`:                     Density of the quadratic matrices
* `random_seed::Int`:                           Random seed used for the Random package
* `is_int_fixed::Bool`:                         Parameter indicating whether all integer variables are fixed (True) or not (False) *
                                                depending on this parameter the correspondent solver is used
* `solver_time_limit::Float64`:                 Time limit for the solver
"""


struct MIP_inital_parameters

    # sceanrios related parameters
    num_scen::Int
    scen_prob::Vector{Float64}

    # parameters per scenario
    num_int_var::Int
    num_cont_var::Int
    num_const::Int

    # parameters for the Random
    quad_mat_dens::Float64
    random_seed::Int

    # parameters for the JuMP
    is_int_fixed::Bool
    solver_time_limit::Float64

end
