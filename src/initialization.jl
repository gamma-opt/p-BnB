
using JuMP, Gurobi, Random, LinearAlgebra, SparseArrays, Suppressor, Ipopt, Profile

# types
include("types/gurobi_solver_parameters.jl")
include("types/bundle_method_parameters.jl")
include("types/dynamic_precision_based_RNMDT_parameters.jl")
include("types/MIP_initial_parameters.jl")
include("types/MIP_generated_parameters.jl")
include("types/node.jl")
include("types/bnb_model.jl")

# utils
include("utils/branching.jl")
include("utils/delta_computation.jl")
include("utils/models_generation.jl")
include("utils/node_generation.jl")
include("utils/node_selection.jl")
include("utils/parameters_generation.jl")

# schemes
include("schemes/bundle_method.jl")
include("schemes/bnb_solve.jl")
include("schemes/dynamic_precision_based_RNMDT_LR_algorithm.jl")

## defifing parameters for gurobi optimizer

Gurobi_Method = 4
Gurobi_OutputFlag = 0
Gurobi_MIPgap = 0
Gurobi_Threads = 1
Gurobi_NonConvex = 2

# generating the structure containing all the parameters values for the Gurobi solver
gurobi_parameters = gurobi_solver_parameters(Gurobi_Method,
                                            Gurobi_OutputFlag,
                                            Gurobi_MIPgap,
                                            Gurobi_Threads,
                                            Gurobi_NonConvex)

## defining general intial parameters for the optimisation problem

g_number_of_scenarios = 4
g_probabilities_of_the_scenarios = [0.10, 0.40, 0.40, 0.10]

# parameters defined per scenario
g_frist_stage_variables_number = 5
g_first_stage_integer_variables_number = 4
g_second_stage_variables_number = 5
g_second_stage_integer_variables_number = 5
g_constraints_number = 2

# parameters for the Random
g_quadratic_matrices_density  = 0.8
g_random_seed = 1

# parameters for the JuMP
g_solver_time_limit = 100 # in seconds

# RNMDT related parameters
g_RNMDT_is_used = true
g_RNMDT_precision_factor = -1 .* ones(g_second_stage_variables_number, g_number_of_scenarios)

## defining parameters for the bundle method

# general parameters
bm_parallelisation_is_used = false

# stoping criteria parameters
bm_max_number_of_iterations = 100
bm_number_of_iteration_for_checking = 4
bm_eps_stop = 0.1

# algorithm performace parameters
bm_m = 0.7
bm_d = 10

# defining  max/min values for the centre of gravity
bm_initial_centre_of_gravity_min = 0
bm_initial_centre_of_gravity_max = 0

# generating the structure containing all the parameters values for bundle method
bm_parameters = bm_input(bm_parallelisation_is_used,
                         bm_max_number_of_iterations,
                         bm_number_of_iteration_for_checking,
                         bm_eps_stop,
                         bm_m,
                         bm_d,
                         bm_initial_centre_of_gravity_min,
                         bm_initial_centre_of_gravity_max)

## defining the parameters for dynamic precision-based RNMDT algorithm


dp_rnmdt_N1_percentage = 10 # the perecnatuge that is used to calculate N1 as some percent of the total number of second stage decision variables in the primal problem
dp_rnmdt_N2 = 5
dp_rnmdt_tolerance = 0.001
dp_rnmdt_time_limit = 200 # in seconds
dp_rndmt_max_number_of_iterations = 100

# generating the structure containing all the parameter values for dynamic precision-based RNMDT algorithm
dp_rnmdt_parameters = dp_based_RNMDT_input(dp_rnmdt_N1_percentage,
                                           dp_rnmdt_N2,
                                           dp_rnmdt_tolerance,
                                           dp_rnmdt_time_limit,
                                           dp_rndmt_max_number_of_iterations)

## generating the structure containing initial parameters

initial_parameters = MIP_initial_parameters( g_number_of_scenarios,
                                             g_probabilities_of_the_scenarios,
                                             g_frist_stage_variables_number,
                                             g_first_stage_integer_variables_number,
                                             g_second_stage_variables_number,
                                             g_second_stage_integer_variables_number,
                                             g_constraints_number,
                                             g_quadratic_matrices_density,
                                             g_random_seed,
                                             g_solver_time_limit,
                                             g_RNMDT_is_used,
                                             g_RNMDT_precision_factor,
                                             bm_parameters,
                                             gurobi_parameters)

## generating the structure containing the constraints and objective related parameters
generated_parameters = parameters_generation(initial_parameters)

## defining the parameters values for BnB algorithm

# non-anticipativity condition - related tolerance
non_ant_tol = 1E-3

# parameter forcing disjunction of the resulting node subproblem feasible sets
# when the resulting node subproblem integrality constraints are violated
tol_bb = 0.1

# integrality condition - related tolerance
integrality_tolerance = 1E-3
