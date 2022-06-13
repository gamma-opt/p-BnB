primal_test_problem = MIP_generation(initial_parameters, generated_parameters)

#@objective(primal_test_problem, Min, -27*primal_test_problem[:y][1,1]^2 - 9*primal_test_problem[:y][2,1]^2 - 35*primal_test_problem[:y][1,2]^2 - 7*primal_test_problem[:y][2,2]^2 - 6.78*primal_test_problem[:x][1,1]
#    - 4.41*primal_test_problem[:x][2,1] - 6.09*primal_test_problem[:y][1,1] - 0.72*primal_test_problem[:y][2,1] - 15.82*primal_test_problem[:x][1,2] - 10.29*primal_test_problem[:x][2,2] - 44.8*primal_test_problem[:y][1,2] - 28.56*primal_test_problem[:y][2,2])

status = optimize!(primal_test_problem)

print(primal_test_problem)

# plot the results
print("\n********************* PRIMAL PROBLEM ********************\n\n")
print("objective function value : $(objective_value(primal_test_problem))\n")
print("first stage variables : $(value.(primal_test_problem[:x]))\n")
print("second stage variables : $(value.(primal_test_problem[:y]))\n")
print("dual gap : $(MOI.get(primal_test_problem, MOI.RelativeGap()))\n")
print("\n**********************************************************\n")

center_of_gravity_min = 0.0
center_of_gravity_max = 0.0


initial_centre_of_gravity = intial_centre_of_gravity_generation(center_of_gravity_min, center_of_gravity_max, initial_parameters.random_seed, initial_parameters.num_scen, initial_parameters.num_first_stage_var)
RNMDT_subproblems = primal_problem_based_lagrangian_relaxation_generation(initial_parameters, generated_parameters)
try_node = node(primal_problem, RNMDT_subproblems , initial_parameters, generated_parameters)

@constraint(try_node.dual_subproblems[1], try_node.dual_subproblems[1][:x][2]>=68.68755662000001)
@constraint(try_node.dual_subproblems[2], try_node.dual_subproblems[2][:x][2]>=68.68755662000001)
@constraint(try_node.dual_subproblems[3], try_node.dual_subproblems[3][:x][2]>=68.68755662000001)
print(try_node.dual_subproblems[1])

result = proximal_bundle_method(try_node, initial_centre_of_gravity, 1e-9)

@show result.dual_objective_value[end]
@show result.variables_values[1]
current_avg_first_stage_var = round.(sum(try_node.initial_parameters.scen_prob[j] .* result.variables_values[1][:,j] for j = 1:try_node.initial_parameters.num_scen), digits = 8)

init_time = time()
output = bnb_solve(initial_parameters, non_ant_tol, tol_bb, integrality_tolerance)
final_time = time()-init_time
print("\n****************** Caroe and Schultz BnB ****************\n\n")
print("objective function value lower bound : $(output.LBDg)\n")
print("objective function value upper bound : $(output.UBDg)\n")
print("first stage variables : $(output.soln_val)\n")
print("dual gap : $((output.UBDg - output.LBDg) / output.LBDg)\n")
print("\n*********************************************************\n")
