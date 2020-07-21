"""
    BnBModel
Stores attributes of stack used to solve BnB problem. Has the following fields:
* `Nodes::Vector{Model}`:                       Nodes associated models stack
* `id::Vector{Int64}`:                          Node ID for each stack item
* `UBDg::Float64`:                              Global Upper Bound
* `LBDg::Float64`:                              Global Lower Bound
* `UBDg_hist::Vector{Float64}`:                 Value history UBD problem
* `LBDg_hist::Vector{Float64}`:                 Value history LBD problem
* `UBDgtime::Vector{Float64}`:                  Run time history UBD problem
* `LBDgtime::Vector{Float64}`:                  Run time history LBD problem
* `max_id::Int64`:                              Max node used
* `soln_hist::Vector{Vector{Float64}}`:         Solution value history
* `soln_val::Float64`:                          Solution value found
* `first_fnd::Bool`:                            Has a solution been found
"""
mutable struct bnb_model
  nodes::Vector{node} # nodes associated models stack
  id::Vector{Int64}  # Node ID
  UBDg::Float64 # Global Upper Bound
  LBDg::Float64 # Global Lower Bound
  UBDg_hist::Vector{Float64} # Value history UBD problem
  LBDg_hist::Vector{Float64} # Value history LBD problem
  UBDg_time_hist::Vector{Float64} # Run time history UBD problem
  LBDg_time_hist::Vector{Float64} # Run time history LBD problem
  max_id::Int64  # Max node used
  soln_hist::Vector{Vector{Float64}}
  soln_val::Vector{Float64}
  first_fnd::Bool
end

"""
    BnBModel(X::Model)
Initializes a `BnBModel` with `.Nodes` = [X with relaxed integrality constraints].
"""
bnb_model(initial_parameters::MIP_initial_parameters) = bnb_model(
                                    [zero_node_generation(initial_parameters)],
                                    [1],
                                    Inf,
                                    -Inf,
                                    [Inf],
                                    [-Inf],
                                    [0.0],
                                    [0.0],
                                    1,
                                    [],
                                    [],
                                    false)
