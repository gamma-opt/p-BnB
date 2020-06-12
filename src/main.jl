using JuMP, Gurobi, Random, LinearAlgebra, SparseArrays

include("types/MIP_generated_parameters.jl")
include("types/MIP_initial_parameters.jl")
include("utils/parameters_generation.jl")
include("utils/models_generation.jl")

initial_parameters = MIP_initial_parameters( 4, [0.25, 0.25, 0.25, 0.25], 3, 3, 1, 0.8, 1, false, 7200 )
model = MIP_generation(initial_parameters)



optimize!(model)

objective_value(model)


l_models = MIP_lagrangian_relaxation_generation(initial_parameters)
optimize!(l_models[1])
objective_value(l_models[1])
optimize!(l_models[2])
objective_value(l_models[2])
optimize!(l_models[3])
objective_value(l_models[3])
optimize!(l_models[4])
objective_value(l_models[4])

objective_value(model) - (objective_value(l_models[1])+objective_value(l_models[2]) + objective_value(l_models[3])+objective_value(l_models[4]))
