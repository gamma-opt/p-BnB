#src_link  =  "/scratch/work/belyakn1/BnB_p_lagrangian/src/"
#src_link  =  "/Users/nikitabelyak/Dropbox (Aalto)/branch-and-bound-caroe-and-schultz/src/"
using JuMP, Gurobi, Random, LinearAlgebra, SparseArrays, Suppressor, Ipopt, Profile, Base.Threads
using DataFrames, XLSX, Dates, Plots

# types
include(src_link*"types/gurobi_solver_parameters.jl")
include(src_link*"types/bundle_method_parameters.jl")
include(src_link*"types/dynamic_precision_based_RNMDT_parameters.jl")
include(src_link*"types/PH_SDM_initial_parameters.jl")
include(src_link*"types/MIP_initial_parameters.jl")
include(src_link*"types/MIP_generated_parameters.jl")
include(src_link*"types/node.jl")
include(src_link*"types/bnb_model.jl")


# utils
include(src_link*"utils/branching.jl")
include(src_link*"utils/delta_computation.jl")
include(src_link*"utils/models_generation.jl")
include(src_link*"utils/node_generation.jl")
include(src_link*"utils/node_selection.jl")
include(src_link*"utils/parameters_generation.jl")
include(src_link*"utils/initialisation_function.jl")
include(src_link*"utils/FW_PH_V_0_initialisation.jl")
include(src_link*"utils/nodes_fathom.jl")
include(src_link*"utils/RNMDT_gap_computation.jl")
include(src_link*"utils/penalty_parameter_update.jl")
include(src_link*"utils/auxiliary_find_repeated_elements_in_V.jl")

# schemes
include(src_link*"schemes/bundle_method.jl")
include(src_link*"schemes/proximal_bundle_method.jl")
include(src_link*"schemes/bnb_solve.jl")
include(src_link*"schemes/dynamic_precision_based_RNMDT_LR_algorithm.jl")
include(src_link*"schemes/FW-PH.jl")
include(src_link*"schemes/simplical_decomposition_method.jl")

## defining the parameters values for BnB algorithm

# non-anticipativity condition - related tolerance
non_ant_tol = 1E-6

# parameter forcing disjunction of the resulting node subproblem feasible sets
# when the resulting node subproblem non-anticipativity constraints are violated
tol_bb = 1E-6

# integrality condition - related tolerance
integrality_tolerance = 1E-6
