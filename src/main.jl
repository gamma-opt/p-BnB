using JuMP, Gurobi, Random, LinearAlgebra, SparseArrays, Suppressor

include("types/MIP_generated_parameters.jl")
include("types/MIP_initial_parameters.jl")
include("types/node.jl")
include("types/bnb_model.jl")
include("utils/parameters_generation.jl")
include("utils/models_generation.jl")

initial_parameters = MIP_initial_parameters( 4, [0.25, 0.25, 0.25, 0.25], 3, 3, 1, 0.8, 1, false, 7200 )

bnb = bnb_model(initial_parameters)
child = child_node_generation(bnb.nodes[1], 1, "<=", 4)
optimize!(child.dual_subproblems[1])


bnb.nodes[1].primal_problem[:x]

model = MIP_generation(initial_parameters)
JuMP.unset_integer.(model[:x])
is_integer.(model[:x])
optimize!(model)
objective_value(model)
value.(model[:x])

l_models = MIP_lagrangian_relaxation_generation(initial_parameters)
optimize!(l_models[1])
objective_value(l_models[1])
x1
optimize!(l_models[2])
objective_value(l_models[2])
optimize!(l_models[3])
objective_value(l_models[3])
optimize!(l_models[4])
objective_value(l_models[4])

objective_value(model) - (objective_value(l_models[1])+objective_value(l_models[2]) + objective_value(l_models[3])+objective_value(l_models[4]))
