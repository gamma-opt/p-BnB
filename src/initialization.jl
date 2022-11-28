
using JuMP, Gurobi, Random, LinearAlgebra, SparseArrays, Suppressor, Ipopt, Profile, Base.Threads
using DataFrames, XLSX, Dates, CSV, MathOptInterface

# types
include(src_link*"types/gurobi_solver_parameters.jl")
include(src_link*"types/bundle_method_parameters.jl")
include(src_link*"types/dynamic_precision_based_RNMDT_parameters.jl")
include(src_link*"types/PH_SDM_initial_parameters.jl")
include(src_link*"types/pooling_problem.jl")
include(src_link*"types/MIP_initial_parameters.jl")
include(src_link*"types/MIP_generated_parameters.jl")
include(src_link*"types/node.jl")
include(src_link*"types/bnb_model.jl")



# utils
include(src_link*"utils/branching.jl")
include(src_link*"utils/delta_computation.jl")
if pp_min_max
    include(src_link*"utils/models_generation.jl")
else
    include(src_link*"utils/models_generation_old.jl")
end
include(src_link*"utils/node_generation.jl")
include(src_link*"utils/node_selection.jl")
if pp_min_max
    include(src_link*"utils/parameters_generation.jl")
else
    include(src_link*"utils/parameters_generation_old.jl")
end
include(src_link*"utils/initialisation_function.jl")
#include(src_link*"utils/initialisation_function_branching.jl")

if pp_min_max
    include(src_link*"utils/FW_PH_V_0_initialisation.jl")
else
    include(src_link*"utils/FW_PH_V_0_initialisation_old.jl")
end
include(src_link*"utils/nodes_fathom.jl")
include(src_link*"utils/RNMDT_gap_computation.jl")
include(src_link*"utils/penalty_parameter_update.jl")
include(src_link*"utils/auxiliary_find_repeated_elements_in_V.jl")
include.(filter(contains(r".jl$"), readdir(src_link*"utils/pooling_problem"; join=true)))
include(src_link*"utils/optimisation_results_print_out.jl")


# schemes
include(src_link*"schemes/bundle_method.jl")
if pp_min_max
    include(src_link*"schemes/proximal_bundle_method.jl")
else
    include(src_link*"schemes/proximal_bundle_method_old.jl")
end
include(src_link*"schemes/bnb_solve.jl")
include(src_link*"schemes/dynamic_precision_based_RNMDT_LR_algorithm.jl")
if pp_min_max
    include(src_link*"schemes/FW-PH.jl")
    include(src_link*"schemes/simplical_decomposition_method.jl")
else
    include(src_link*"schemes/FW-PH_old.jl")
    include(src_link*"schemes/simplical_decomposition_method_old.jl")
end
#


## defining the parameters values for BnB algorithm

# non-anticipativity condition - related tolerance
non_ant_tol = 1E-6

# parameter forcing disjunction of the resulting node subproblem feasible sets
# when the resulting node subproblem non-anticipativity constraints are violated
tol_bb = 1E-6

# integrality condition - related tolerance
integrality_tolerance = 1E-6
