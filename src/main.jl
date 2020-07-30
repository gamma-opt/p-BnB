

using JuMP, Gurobi, Random, LinearAlgebra, SparseArrays, Suppressor, Ipopt

include("types/MIP_initial_parameters.jl")
include("types/MIP_generated_parameters.jl")
include("types/node.jl")
include("types/bnb_model.jl")
include("types/bundle_method_output.jl")

include("utils/branching.jl")
include("utils/delta_computation.jl")
include("utils/models_generation.jl")
include("utils/node_generation.jl")
include("utils/node_selection.jl")
include("utils/parameters_generation.jl")
include("utils/update_LB.jl")

include("schemes/bundle_method.jl")
include("schemes/bnb_solve.jl")

## defining intial parameters for the optimisation problem
initial_parameters = MIP_initial_parameters( 4, [0.15, 0.45, 0.35, 0.05], 4, 4, 4, 0.8, 1, false, 100 )

## solving the full-scale problem using Gurobi solver with predefined time limit
model = MIP_generation(initial_parameters)
optimize!(model)
objective_value(model)
value.(model[:x])
gap = MOI.get(model, MOI.RelativeGap())

## solving the problem using caroe and schultz bnb method with predefined parameters
ouput = bnb_solve(initial_parameters, 0.1, 0.1)
(output.UBDg - output.LBDg) /  output.LBDg
