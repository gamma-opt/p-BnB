using JuMP, Gurobi, Random, LinearAlgebra, SparseArrays

include("types/MIP_generated_parameters.jl")
include("types/MIP_initial_parameters.jl")
include("utils/parameters_generation.jl")
include("utils/models_generation.jl")

inital_parameters = MIP_inital_parameters( 2, [0.5, 0.5], 3, 3, 1, 0.8, 1, false, 7200 )
