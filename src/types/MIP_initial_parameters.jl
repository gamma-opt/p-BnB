"""
    MIP_initial_parameters
Stores attributes for generating MIP JuMP model. Has the following fields:
* `num_scen::Int`:                                  Number of scenarios used in the model
* `scen_prob::Vector{Float64}`:                     Vector of probabilities associated with each scenario
* `num_fist_stage_var::Int`:                        Number of first stage variables associated with each scenario
*` num_int_var_first_stage`:                        Number of integer variables among first stage varaibles associated with each scenario
* `num_second_stage_var::Int`:                      Number of second stage variables associated with each scenario
*` num_int_var_second_stage`:                       Number of integer variables among second stage varaibles associated with each scenario
* `num_const::Int`:                                 Number of constraints associated with each scenario
* `quad_mat_dens::Float64`:                         Density of the quadratic matrices
* `random_seed::Int`:                               Random seed used for the Random package

* `μ::Float64`:                                     Big enough constant use as a penalty parameter for the slack variables (to ensure full recourse)

* `al_is_used::Bool`:                               Parameter that takes value TRUE if Augemnted lagrangian is used and FALSE otherwise
* `al_penalty_parameter::Float64`:                  Augemnted lagrangian relaxation penalty parameter

* `solver_time_limit::Float64`:                     Time limit for the solver

* `RNMDT_is_used::Bool`:                            Is RNMDT technique used to appriximate MIQCQP with MIP
* `RNMDT_precision_factor:: Array{Int}`:            RNMDT-based psreicsion factor values for the second stage variables


* `bm_parameters::bm_input`:                        Parameters for bundle method

* `PH_SDM_parameters::PH_SDM_input`:                Parameters for Frank-Wolfe Progressive Hedging and Simplical Decomposition Method

* `gurobi_parameters::gurobi_solver_parameters`:    Parameters for gurobi optimizer
"""


mutable struct MIP_initial_parameters

    # sceanrios related parameters
    num_scen::Int
    scen_prob::Vector{Float64}

    # parameters per scenario
    num_first_stage_var::Int
    num_int_var_first_stage::Int
    num_second_stage_var::Int
    num_int_var_second_stage::Int
    num_const::Int

    # parameters for the Random
    quad_mat_dens::Float64
    random_seed::Int

    # slack variables related parameters (ensuring full recourse)
    μ::Float64

    # augmented lagrangian related parameters
    al_is_used::Bool
    al_penalty_parameter::Float64

    # parameters for the JuMP
    solver_time_limit::Float64

    # parameters for RNMDT
    RNMDT_is_used::Bool
    RNMDT_precision_factor:: Array{Int}

    # parameters for bundle method
    bm_parameters::bm_input

    # Frank-Wolfe Progressive Hedging and Simplical Decomposition Method related parameters
    PH_SDM_parameters::PH_SDM_input

    # Gurobi related parameters
    gurobi_parameters::gurobi_solver_parameters

end
