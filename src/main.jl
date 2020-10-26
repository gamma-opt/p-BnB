include("initialization.jl")

## solving the full-scale problem using Gurobi solver with predefined time limit

# generating primal problem
primal_problem = MIP_generation(initial_parameters, generated_parameters)

# fix the values of the first-stage decision variables if needed
values_to_be_fixed = [64, 0, 1, 73, 0]
for s = 1:initial_parameters.num_scen
    for j = 1:initial_parameters.num_int_var_first_stage
        fix(primal_problem[:x][j,s], values_to_be_fixed[j])
    end
end

# optimising the problem
optimize!(primal_problem)

# plot the results
print("\n********************* PRIMAL PROBLEM ********************\n\n")
print("objective function value : $(objective_value(primal_problem))\n")
print("first stage variables : $(value.(primal_problem[:x]))\n")
print("second stage variables : $(value.(primal_problem[:y]))\n")
print("dual gap : $(MOI.get(primal_problem, MOI.RelativeGap()))\n")
print("\n**********************************************************\n")





## solving RNMDT based relaxation with fixed precision factor

# generating RNMDT relxation
RNMDT_relaxation = RNMDT_based_problem_generation(initial_parameters, generated_parameters)

# optimising the problem
optimize!(RNMDT_relaxation)

# plot the results
print("\n******************** RNMDT relaxation *******************\n\n")
print("objective function value : $(objective_value(RNMDT_relaxation))\n")
print("first stage variables : $(value.(RNMDT_relaxation[:x]))\n")
print("second stage variables : $(value.(RNMDT_relaxation[:y]))\n")
print("dual gap : $(MOI.get(RNMDT_relaxation, MOI.RelativeGap()))\n")
print("\n*********************************************************\n")

## solving the problem using caroe and schultz bnb method with predefined parameters
init_time = time()
output = bnb_solve(initial_parameters, non_ant_tol, tol_bb, integrality_tolerance)
final_time = time()-init_time

# plot the results
print("\n****************** Caroe and Schultz BnB ****************\n\n")
print("objective function value lower bound : $(output.LBDg)\n")
print("objective function value upper bound : $(output.UBDg)\n")
print("first stage variables : $(output.soln_val)\n")
print("dual gap : $((output.UBDg - output.LBDg) / output.LBDg)\n")
print("\n*********************************************************\n")

## checking the precision of the bundle method

# generating the inital values for the center of gravity
center_of_gravity_min = 0
center_of_gravity_max = 0

bnb_struct = bnb_model(initial_parameters)
dynamic_precision_RNMDT_algorithm(bnb_struct.nodes[1])

s = primal_problem_based_lagrangian_relaxation_generation(initial_parameters, generated_parameters)

print(s[1])

initial_centre_of_gravity = intial_centre_of_gravity_generation(initial_parameters.bm_parameters.initial_centre_of_gravity_min, initial_parameters.bm_parameters.initial_centre_of_gravity_max, initial_parameters.random_seed, initial_parameters.num_scen, initial_parameters.num_first_stage_var)
RNMDT_subproblems = RNMDT_based_lagrangian_relaxation_problem_generation(initial_parameters, generated_parameters, initial_parameters.RNMDT_precision_factor)
try_node = node(primal_problem, RNMDT_subproblems , initial_parameters, generated_parameters)
bundle_method(try_node, initial_centre_of_gravity)
