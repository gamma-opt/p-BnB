function proximal_bundle_method( vector_of_subproblems::Vector{Model}, initial_penality_parameter::Float64)



k = 0  # iteration counter
lagrangian_multipliers = Vector{Vector{Float64}}

## zero iteration

# bundle_subproblem_1
m = Vector{Model(Gurobi.Optimizer)}
@variable(m, theta[1 : initial_parameters.num_scen])
@objective(m, Max, sum(theta))
