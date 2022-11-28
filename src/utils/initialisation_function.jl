
function initialisation(g_number_of_scenarios, g_frist_stage_variables_number, g_second_stage_variables_number, g_constraints_number, precision_p, AL, r_seed )


## defifing parameters for gurobi optimizer

Gurobi_Method = -1
Gurobi_OutputFlag = 0
#Gurobi_MIPgap = 0
Gurobi_IntFeasTol = 1E-6
Gurobi_FeasibilityTol = 1E-6
Gurobi_OptimalityTol = 1E-6
Gurobi_Threads = 1
Gurobi_NonConvex = 2
Numeric_Focus = 0

# generating the structure containing all the parameters values for the Gurobi solver
gurobi_parameters = gurobi_solver_parameters(Gurobi_Method,
                                            Gurobi_OutputFlag,
                                            #Gurobi_MIPgap,
                                            Gurobi_IntFeasTol,
                                            Gurobi_FeasibilityTol,
                                            Gurobi_OptimalityTol,
                                            Gurobi_Threads,
                                            Gurobi_NonConvex,
                                            Numeric_Focus)

## defining general intial parameters for the optimisation problem

# parameters for the Random
g_quadratic_matrices_density  = 0.9
g_random_seed = r_seed
# defining max and min values for the quadratic matirces coefficients in objective
g_obj_quad_min = -100
g_obj_quad_max = 100
# defining max and min values for first-stage variables linear coefficients in objective
g_obj_lin_fs_min = 0
g_obj_lin_fs_smax = 100
# defining max and min values for second-stage variables linear coefficients in objective
g_obj_lin_ss_min = -100 
g_obj_lin_ss_smax = 100
# defining max and min values for the quadratic matirces coefficients in constraints
g_con_quad_min = 0
g_con_quad_max = 1000
# defining max and min values for first-stage variables linear coefficients in constraints
g_con_lin_fs_min = 0
g_con_lin_fs_max = 100
# defining max and min values for second-stage variables linear coefficients in constraints
g_con_lin_ss_min = -100
g_con_lin_ss_max = 100
# defining max and min value for the affine constant in the constraints
g_con_aff_min = 1000 
g_con_aff_max = 100000
# defining whether the first-stage veraible are binary (True) or not (False)
g_bin_con_fs = false
# defining max and min values for the first-stage variables' box constraints
g_box_con_fs_min = 0 
g_box_con_fs_max = 10
# defining max and min values for the second-stage variables' box constraints
g_box_con_ss_min = 0
g_box_con_ss_max = 10



#g_number_of_scenarios = 5
#g_probabilities_of_the_scenarios = [0.05, 0.1, 0.05, 0.2, 0.1, 0.05, 0.1, 0.05, 0.2, 0.1]
Random.seed!(g_random_seed)
g_probabilities_of_the_scenarios = rand(g_number_of_scenarios)
g_probabilities_of_the_scenarios = g_probabilities_of_the_scenarios./sum(g_probabilities_of_the_scenarios)

# parameters defined per scenario
#g_frist_stage_variables_number = 5
g_first_stage_integer_variables_number = g_frist_stage_variables_number
#g_second_stage_variables_number = 5
g_second_stage_integer_variables_number = 0

# slack variables realted parameters
g_μ = 10000000 # penality pramemeter for the objective

#-----------------------------------------------------------------------------
# pooling problem parameters
g_pool_problem_is_used = false

pp_primal_pool_problem_link = chop(src_link, tail = 4) * "pooling_problem_data"
pp_print_out_variables_correlation_map = true
pp_start_value_pooling_nodes = [1.0]
g_flow_cost = 10
g_pool_cost = 10

# stochastic pooling problem parameters
stoch_pp_pool_max_capacity = 10.0
stoch_pp_flow_max_capacity = 10.0
stoch_pp_demand_min_value = 0.0 
stoch_pp_demand_max_value = 5.0
stoch_pp_source_node_min_value = 1.0
stoch_pp_source_node_max_value = 3.0

g_stoc_ppm = stochastic_pooling_problem_parameters(stoch_pp_pool_max_capacity, stoch_pp_flow_max_capacity, stoch_pp_demand_min_value, stoch_pp_demand_max_value, stoch_pp_source_node_min_value, stoch_pp_source_node_max_value)

g_pooling_problem_parameters = pooling_problem_parameters(pp_primal_pool_problem_link, pp_print_out_variables_correlation_map, pp_start_value_pooling_nodes, g_flow_cost, g_pool_cost, g_stoc_ppm)
#-----------------------------------------------------------------------------

# augmented lagrangian parameters
g_al_is_used = AL
# penalty parameter related parameters
g_al_is_fixed = false
g_al_start_value = 100.0
g_al_penalty_parameter = g_al_start_value * ones(g_frist_stage_variables_number)

# parameters for the JuMP
g_solver_time_limit = 3600 # in seconds

# RNMDT related parameters
g_RNMDT_is_used = true
g_RNMDT_precision_factor = precision_p .* ones(g_second_stage_variables_number, g_number_of_scenarios)

## defining parameters for the bundle method

# general parameters
bm_parallelisation_is_used = true

# stoping criteria parameters
bm_max_number_of_iterations = 1000
bm_number_of_iteration_for_checking = 4
bm_eps_stop = 1E-3

# algorithm performace parameters
bm_mR = 0.7 # in [0.5, 1) - stepsize related
bm_mL = 0.3 # in (0, 0.5) - serious step related
bm_min_ssv = 0.001
bm_max_ssv = 1000.0

# defining  max/min values for the centre of gravity
bm_initial_centre_of_gravity_min = 0.0
bm_initial_centre_of_gravity_max = 0.0

# generating the structure containing all the parameters values for bundle method
bm_parameters = bm_input(bm_parallelisation_is_used,
                         bm_max_number_of_iterations,
                         bm_number_of_iteration_for_checking,
                         bm_eps_stop,
                         bm_mR,
                         bm_mL,
                         bm_min_ssv,
                         bm_max_ssv,
                         bm_initial_centre_of_gravity_min,
                         bm_initial_centre_of_gravity_max)
## defining parameters for Frank-Wolfe Progressive Hedging (FW-PH) and Simplical Decomposition Method (SDM)

# defining parameter that affects initial linearisation point of the SDM method
FW_PH_α = 1.0

# defining the tolerance used for stopping criterion in FW-PH
FW_PH_tolerance = 1E-3

# defining the maximum number of iterations used in FW-PH
PH_max_iter = 1000

# The indicator whether we should or not plot the FW-PH porgrees for the root node of the BnB
FW_PH_plot = false

# the maximum time allowed for FW-PH to solve one instance
FW_PH_PH_max_time = 3600

# defining the tolerance used for stopping criterion in SDM
SDM_tolerance = 0.0

# defining the maximum number of iterations used in SDM
SDM_max_iter = 1

PH_SDM_parameters = PH_SDM_input(
                        FW_PH_α,
                        FW_PH_tolerance,
                        PH_max_iter,
                        FW_PH_plot,
                        FW_PH_PH_max_time,
                        SDM_tolerance,
                        SDM_max_iter
                        )

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
                                             g_obj_quad_min,
                                             g_obj_quad_max,
                                             g_obj_lin_fs_min,
                                             g_obj_lin_fs_smax,
                                             g_obj_lin_ss_min,
                                             g_obj_lin_ss_smax,
                                             g_con_quad_min,
                                             g_con_quad_max,
                                             g_con_lin_fs_min,
                                             g_con_lin_fs_max,
                                             g_con_lin_ss_min,
                                             g_con_lin_ss_max,
                                             g_con_aff_min,
                                             g_con_aff_max,
                                             g_bin_con_fs,
                                             g_box_con_fs_min ,
                                             g_box_con_fs_max,
                                             g_box_con_ss_min,
                                             g_box_con_ss_max, 
                                             g_μ,
                                             g_pool_problem_is_used,
                                             g_pooling_problem_parameters,
                                             g_al_is_used,
                                             g_al_is_fixed,
                                             g_al_start_value,
                                             g_al_penalty_parameter,
                                             g_solver_time_limit,
                                             g_RNMDT_is_used,
                                             g_RNMDT_precision_factor,
                                             bm_parameters,
                                             PH_SDM_parameters,
                                             gurobi_parameters)

## generating the structure containing the constraints and objective related parameters
#generated_parameters = parameters_generation(initial_parameters)

    return initial_parameters

end
