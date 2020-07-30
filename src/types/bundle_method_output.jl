"""
    bundle_method_output
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

struct bundle_method_output
    dual_objective_value::Vector{Float64}
    int_var:: Array{Float64}
    cont_var:: Array{Float64}
end
