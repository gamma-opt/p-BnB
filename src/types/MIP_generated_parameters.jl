"""
    MIP_initial_parameters
Stores generated parameters for the MIP JuMP model. Has the following fields:
* `constraint_Qs::Array{Any}`:                  Vector of quadratic matrices for the left hand side of the contraints associated with each scenario
* `constraint_fs::Array{Any}`:                  Vector of affine functions' coefficients for the left hand side of the constraints associated with each scenario
* `objective_Qs::Array{Any}`:                   Vector of quadratic matrices associated with each scenario for the objective function
* `objective_fs::Array{Any}`:                   Vector of affine functions' coefficients for continuous variables associated with each scenario for the objective function
* `objective_c::Array{Float64}`:                Vector of affine functions' coefficients for integer variables for the objective function
* `x_boundaries::Array{Int}`:               Array of box constraints coefficients for integer variables
* `y_boundaries::Array{Float64}`:               Array of box constraints coefficients for continuous variables
"""


struct MIP_generated_parameters

    constraint_Qs::Array{Any}
    constraint_fs::Array{Any}
    objective_Qs::Array{Any}
    objective_fs::Array{Any}
    objective_c::Array{Float64}
    x_boundaries::Array{Int}
    y_boundaries::Array{Float64}

end
