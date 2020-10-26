"""
    MIP_initial_parameters
Stores generated parameters for the MIP JuMP model. Has the following fields:
* `constraint_Qs::Array{Any}`:                  Vector of quadratic matrices for the left hand side of the contraints associated with each scenario
* `constraint_fs::Array{Any}`:                  Vector of affine functions' coefficients for the left hand side of the constraints associated with each scenario
* `objective_Qs::Array{Any}`:                   Vector of quadratic matrices associated with each scenario for the objective function
* `objective_fs::Array{Any}`:                   Vector of affine functions' coefficients for continuous variables associated with each scenario for the objective function
* `objective_c::Array{Float64}`:                Vector of affine functions' coefficients for integer variables for the objective function
* `x_boundaries::Array{Int}`:                   Array of box constraints coefficients for first stage variables (for all scenarios)
* `x_int_indexes::Array{Int}`:                  Array of integer indexes for first stage variables (for all scenarios)
* `x_cont_indexes::Array{Int}`:                 Array of continuous indexes for first stage variables (for all scenarios)
* `y_boundaries::Array{Float64}`:               Array of box constraints coefficients for second stage variables (for each scenario)
* `y_int_indexes::Array{Int}`:                  Array of integer indexes for second stage variables (for each scenario)
* `y_cont_indexes::Array{Int}`:                 Array of continuous indexes for second stage variables (for each scenario)

"""


struct MIP_generated_parameters

    constraint_Qs::Array{Any}
    constraint_fs::Array{Any}
    objective_Qs::Array{Any}
    objective_fs::Array{Any}
    objective_c::Array{Float64}
    x_boundaries::Array{Int}
    x_int_indexes::Array{Int}
    x_cont_indexes::Array{Int}
    y_boundaries::Array{Any}
    y_int_indexes::Array{Int,2}
    y_cont_indexes::Array{Int,2}


end
