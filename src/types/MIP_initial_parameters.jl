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
* `obj_quad_min::Float64`:                          Min value for the quadratic matirces coefficients in objective
* `obj_quad_max::Float64`:                          Max value for the quadratic matirces coefficients in objective
* `obj_lin_fs_min::Float64`:                        Min value for first-stage variables linear coefficients in objective
* `obj_lin_fs_max::Float64`:                        Max value for first-stage variables linear coefficients in objective
* `obj_lin_ss_min::Float64`:                        Min value for second-stage variables linear coefficients in objective
* `obj_lin_ss_smax::Float64`:                       Max value for second-stage variables linear coefficients in objective    
* `con_quad_min::Float64`:                          Min value for the quadratic matirces coefficients in constraints
* `con_quad_max::Float64`:                          Max value for the quadratic matirces coefficients in constraints
* `con_lin_fs_min::Float64`:                        Min value for first-stage variables linear coefficients in constraints
* `con_lin_fs_max::Float64`:                        Max value for first-stage variables linear coefficients in constraints
* `con_lin_ss_min::Float64`:                        Min value for second-stage variables linear coefficients in constraints
* `con_lin_ss_max::Float64`:                        Max value for second-stage variables linear coefficients in constraints     
* `con_aff_min::Float64`:                           Min value for for the affine constant in constraints
* `con_aff_max::Float64`:                           Max value for for the affine constant in constraints         
* `bin_con_fs:: Bool`:                              Indicator whether the first-stage veraible are binary (True) or not (False)       
* `box_con_fs_min::Float64`:                        Min value for the first-stage variables' box constraints     
* `box_con_fs_max::Float64`:                        Max value for the first-stage variables' box constraints       
* `box_con_ss_min::Float64`:                        Min value for the second-stage variables' box constraints       
* `box_con_ss_max::Float64`:                        Max value for the second-stage variables' box constraints         
    
* `μ::Float64`:                                     Big enough constant use as a penalty parameter for the slack variables (to ensure full recourse)

* `pool_problem_is_used::Bool`:                     Is pooling problem used to generate the instances
* `pool_prob_par::pooling_problem_parameters`:      Pooling problem parameters

* `al_is_used::Bool`:                               Parameter that takes value TRUE if Augemnted lagrangian is used and FALSE otherwise
* `al_is_fixed::Bool`:                              The parameter that takes value TRUE in case the Augmented lagrangian penalty term should be fixed to one value (no penalty update procedure allowed)
* `al_start_value::Float64`:                        Augmented lagrangian relaxation penalty parameter starting value that will be spanned to all components        
* `al_penalty_parameter::Array{Float64}`:           Augemnted lagrangian relaxation penalty parameter vector (with the penalty value for each component)

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
    obj_quad_min::Float64
    obj_quad_max::Float64
    obj_lin_fs_min::Float64
    obj_lin_fs_smax::Float64
    obj_lin_ss_min::Float64
    obj_lin_ss_max::Float64
    con_quad_min::Float64
    con_quad_max::Float64
    con_lin_fs_min::Float64
    con_lin_fs_max::Float64
    con_lin_ss_min::Float64
    con_lin_ss_max::Float64
    con_aff_min::Float64
    con_aff_max::Float64
    bin_con_fs:: Bool
    box_con_fs_min::Float64
    box_con_fs_max::Float64
    box_con_ss_min::Float64
    box_con_ss_max::Float64


    # slack variables related parameters (ensuring full recourse)
    μ::Float64

    # pooling problem related parameters
    pool_problem_is_used::Bool
    pool_prob_par::pooling_problem_parameters

    # augmented lagrangian related parameters
    al_is_used::Bool
    al_is_fixed::Bool
    al_start_value::Float64
    al_penalty_parameter::Array{Float64}

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
