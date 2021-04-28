include("initialization.jl")


## solving the full-scale problem using Gurobi solver with predefined time limit

# generating primal problem
primal_problem = MIP_generation(initial_parameters, generated_parameters)

# fix the values of the first-stage decision variables if needed
values_to_be_fixed = [6.0 6.0 6.0; 6.0 6.0 6.0]
for s = 1:initial_parameters.num_scen
    for j = 1:initial_parameters.num_int_var_first_stage
       fix(primal_problem[:x][j,s], values_to_be_fixed[j,s])
    end
end

fix.(primal_problem[:x], x2)

print(primal_problem)
# optimising the problem
@time status = optimize!(primal_problem)

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
fix.(RNMDT_relaxation[:x][:,2],[5, 5])
print(RNMDT_relaxation)

# optimising the problem
init = time()
@time optimize!(RNMDT_relaxation)
finlt = init - time()

# plot the results
print("\n******************** RNMDT relaxation *******************\n\n")
print("objective function value : $(objective_value(RNMDT_relaxation))\n")
print("first stage variables : $(value.(RNMDT_relaxation[:x]))\n")
print("second stage variables : $(value.(RNMDT_relaxation[:y]))\n")
print("dual gap : $(MOI.get(RNMDT_relaxation, MOI.RelativeGap()))\n")
print("\n*********************************************************\n")
print(value.(RNMDT_relaxation[:w_RNMDT]))
10
value.(RNMDT_relaxation[:w_RNMDT])
value.(RNMDT_relaxation[:y])
RNMDT_gap_computation( value.(RNMDT_relaxation[:y]), value.(RNMDT_relaxation[:w_RNMDT]))

## solving primal problem based McCormic relaxation

McCormic_relaxation = primal_problem_based_McCormic_relaxation(initial_parameters, generated_parameters)

# optimising the problem
optimize!(McCormic_relaxation)

# plot the results
print("\n******************** McCormic relaxation *******************\n\n")
print("objective function value : $(objective_value(McCormic_relaxation))\n")
print("first stage variables : $(value.(McCormic_relaxation[:x]))\n")
print("second stage variables : $(value.(McCormic_relaxation[:y]))\n")
print("dual gap : $(MOI.get(McCormic_relaxation, MOI.RelativeGap()))\n")
print("\n*********************************************************\n")

print(McCormic_relaxation)

## solving the problem using caroe and schultz bnb method with predefined parameters
init_time = time()
output = bnb_solve(initial_parameters, non_ant_tol, tol_bb, integrality_tolerance)
final_time = time()-init_time

Juno.profiler()
Profile.clear()

initial_parameters.al_is_used = false

# plot the results
print("\n****************** Caroe and Schultz BnB ****************\n\n")
print("objective function value lower bound : $(output.LBDg)\n")
print("objective function value upper bound : $(output.UBDg)\n")
print("first stage variables : $(output.soln_val)\n")
print("dual gap : $((output.UBDg - output.LBDg) / output.LBDg)\n")
print("\n*********************************************************\n")





## defining the parameters values for BnB algorithm

# non-anticipativity condition - related tolerance
non_ant_tol = 1E-6

# parameter forcing disjunction of the resulting node subproblem feasible sets
# when the resulting node subproblem non-anticipativity constraints are violated
tol_bb = 1E-6

# integrality condition - related tolerance
integrality_tolerance = 1E-6

## Profiling
Profile.clear()
@time bnb_solve(initial_parameters, non_ant_tol, tol_bb, integrality_tolerance)
Juno.profiler(; C = true)

Profile.print()
Profile.print(format=:flat)


@profiler bnb_solve(initial_parameters, non_ant_tol, tol_bb, integrality_tolerance) combine = true

max(output.LBDg_hist)
## checking the precision of the bundle method

# generating the inital values for the center of gravity
center_of_gravity_min = 0.0
center_of_gravity_max = 0.0

initial_centre_of_gravity = intial_centre_of_gravity_generation(center_of_gravity_min, center_of_gravity_max, initial_parameters.random_seed, initial_parameters.num_scen, initial_parameters.num_first_stage_var)
#subproblems = primal_problem_based_lagrangian_relaxation_generation(initial_parameters, generated_parameters)
RNMDT_subproblems = RNMDT_based_lagrangian_relaxation_problem_generation(initial_parameters, generated_parameters, initial_parameters.RNMDT_precision_factor)
try_node = node(primal_problem, RNMDT_subproblems , initial_parameters, generated_parameters)

@constraint(try_node.dual_subproblems[1], try_node.dual_subproblems[1][:x][2]<=3.19726787)
@constraint(try_node.dual_subproblems[2], try_node.dual_subproblems[2][:x][2]<=3.19726787)
@constraint(try_node.dual_subproblems[3], try_node.dual_subproblems[3][:x][2]<=3.19726787)
@constraint(try_node.dual_subproblems[4], try_node.dual_subproblems[4][:x][4]<=3.32056313)
result = proximal_bundle_method(try_node, initial_centre_of_gravity, 1e-9)

print(try_node.dual_subproblems[2])
dual_subproblems_values(try_node, initial_centre_of_gravity)

model = try_node.dual_subproblems[4]
optimize!(model)
Int(termination_status(model))
typeof(sss)
Int(sss)


JuMP.compute_conflict!(model)
computeIIS(model)
grb_model=backend(model)
grb.write()

MOI.get(model, Gurobi.ModelAttribute("LB"))
MOI.get(model.moi_backend, Gurobi.ConstraintConflictStatus(), con1.index)


@show initial_parameters

@show !isempty(generated_parameters.x_int_indexes)

##


McCormic_relaxation = primal_problem_based_McCormic_relaxation(initial_parameters, generated_parameters)
@constraint(McCormic_relaxation,
    4.77695 .<= McCormic_relaxation[:y][1, 1] .<= 6.2303)

optimize!(McCormic_relaxation)
objective_value(McCormic_relaxation)
value.(McCormic_relaxation[:x])

tighten_bounds = bounds_tightening(initial_parameters, generated_parameters, McCormic_relaxation, -13632.502629213033)
 @show tighten_bounds

generated_parameters.y_boundaries = tighten_bounds

RNMDT_relaxation = RNMDT_based_problem_generation(initial_parameters, generated_parameters)
Result_x = Array{Array{Float64}}(undef, 4, 3)
s_initial_parameters = initialsition(9, 9, 9, 9, -2)
s_generated_parameters = parameters_generation(s_initial_parameters)
fix.(RNMDT_relaxation[:x], x1)
RNMDT_relaxation = RNMDT_based_problem_generation(s_initial_parameters, s_generated_parameters)
status = optimize!(RNMDT_relaxation)

f1 = objective_value(RNMDT_relaxation)
f2 = objective_value(RNMDT_relaxation)




    - sum( s_initial_parameters.scen_prob[s] *
        (
        sum(s_generated_parameters.objective_Qs[s][i, j] * abs(w1[i, j, s] - w2[i, j, s])
                for i = 1 : s_initial_parameters.num_second_stage_var,
                    j = 1 : s_initial_parameters.num_second_stage_var)
            + sum( abs(y1[j, s] - y2[j, s]) * s_generated_parameters.objective_fs[s][j] for j = 1:s_initial_parameters.num_second_stage_var)
        )
    for s = 1 : s_initial_parameters.num_scen)
# 26
x1 = value.(RNMDT_relaxation[:x])
x5 = value.(RNMDT_relaxation[:x])
x6 = value.(RNMDT_relaxation[:x])
x7 = value.(RNMDT_relaxation[:x])
x8 = value.(RNMDT_relaxation[:x])
x9 = value.(RNMDT_relaxation[:x])
x10 = value.(RNMDT_relaxation[:x])
x11 = value.(RNMDT_relaxation[:x])
x12 = value.(RNMDT_relaxation[:x])
y1 = value.(RNMDT_relaxation[:y])
w1 = value.(RNMDT_relaxation[:w_RNMDT])
x2 = value.(RNMDT_relaxation[:x])
y2 = value.(RNMDT_relaxation[:y])
w2 = value.(RNMDT_relaxation[:w_RNMDT])

p_12 = y1 .- y2
w_12 = w1 .- w2
w_12[:,:,1]

norm_w_12 = sum(norm(w_12[:, j, s]) for j = 1:3, s = 1:3)
norm_y_12 = sum(norm(p_12[:,s]) for s = 1:3 )

p_23 = z2 .- z1
count_s = 0
for s = [3, 6, 9, 12]
    count_s = count_s+1
    for p = -1:-1:-3
        s_initial_parameters = initialsition(s, s, s, s, p)
        s_generated_parameters = parameters_generation(s_initial_parameters)
        RNMDT_relaxation = RNMDT_based_problem_generation(s_initial_parameters, s_generated_parameters)
        status = optimize!(RNMDT_relaxation)
        Result_x[count_s,-p] = value.(RNMDT_relaxation[:x])
    end
end
status = optimize!(RNMDT_relaxation)

Result_x[1,:]

# plot the results
print("\n******************** RNMDT relaxation *******************\n\n"),0,0
print("objective function value : $(objective_value(RNMDT_relaxation))\n")
print("first stage variables : $(value.(RNMDT_relaxation[:x]))\n")
print("second stage variables : $(value.(RNMDT_relaxation[:y]))\n")
print("dual gap : $(MOI.get(RNMDT_relaxation, MOI.RelativeGap()))\n")
print("\n*********************************************************\n")

w = [value.(RNMDT_relaxation[:w_RNMDT][i,j,s]) .-  value.(RNMDT_relaxation[:y])[i,s]* value.(RNMDT_relaxation[:y])[j,s] for i = 1:initial_parameters.num_second_stage_var, j = 1:initial_parameters.num_second_stage_var, s= 1:initial_parameters.num_scen]
@show w[1,1,2]


# generating the inital values for the center of gravity
center_of_gravity_min = 0.0
center_of_gravity_max = 0.0

initial_centre_of_gravity = intial_centre_of_gravity_generation(center_of_gravity_min, center_of_gravity_max, initial_parameters.random_seed, initial_parameters.num_scen, initial_parameters.num_first_stage_var)
#subproblems = primal_problem_based_lagrangian_relaxation_generation(initial_parameters, generated_parameters)
RNMDT_subproblems = RNMDT_based_lagrangian_relaxation_problem_generation(initial_parameters, generated_parameters, initial_parameters.RNMDT_precision_factor)
try_node = node(primal_problem, RNMDT_subproblems , initial_parameters, generated_parameters)
try_node = zero_node_generation(initial_parameters)

fix.(RNMDT_subproblems[1][:x], [6, 6 , 10])
optimize!(RNMDT_subproblems[1])
value.(RNMDT_subproblems[1][:y])
value.(RNMDT_subproblems[1][:w_RNMDT])
y_0 = value.(RNMDT_relaxation[:y])
w_RNMDT_0 = value.(RNMDT_relaxation[:w_RNMDT])

initial_centre_of_gravity[1]

print(try_node.dual_subproblems[2])
x_SDM = [5.0 5.0]
y_SDM = [6.0 6.258]
w_SDM = [27.0 45.78782112274024; 43.14391056137012 37.57564224548049]
ss = Array{Array{Any}}(undef,2)
ss = [x_SDM, y_SDM, w_SDM]

SDM(2, try_node, ss, x_SDM, [0.0 0.0], [6.0 6.0], 10, 0.0001)

x_0 = value.(RNMDT_relaxation[:x])
y_0 = value.(RNMDT_relaxation[:y])
w_RNMDT_0 = value.(RNMDT_relaxation[:w_RNMDT])
z_0 = value.(RNMDT_relaxation[:z])
w_RNMDT_0[:,:,1] = value.(RNMDT_subproblems[1][:w_RNMDT])
w_RNMDT_0[:,:,2] = value.(RNMDT_subproblems[2][:w_RNMDT])
w_RNMDT_0[:,:,2] = value.(RNMDT_subproblems[2][:w_RNMDT])
V_0 = Array{AbstractArray{Vector{Array{Float64}},1}}(undef, initial_parameters.num_scen)
[ V_0[s] = [[ x_0[:,s], y_0[:,s], w_RNMDT_0[:,:,s], z_0[s,:] ]] for s = 1:initial_parameters.num_scen]

V_0, x_0 = FW_PH_V_0_initialisation( try_node, initial_centre_of_gravity)
V_0[1]
V_0[5]
x_0[2]
x_0
dual_values, x_0, y_0, w_RNMDT_0, z_0 = dual_subproblems_values(try_node, initial_centre_of_gravity)

V_0 = Array{Any}(undef, initial_parameters.num_scen)
[ V_0[s] = [x_0[:,s], y_0[:,s], w_RNMDT_0[:,:,s] ] for s = 1:initial_parameters.num_scen]

x_k, y_k, w_RNMDT_k, z_k, w_k, ϕ_k_bundle = FW_PH_bundle(try_node, V_0, x_0, y_0, w_RNMDT_0, initial_centre_of_gravity)
x_k, y_k, w_RNMDT_k, z_k, w_k, ϕ_k = FW_PH(try_node, V_0, x_0, initial_centre_of_gravity)
ss = [x_k[:]
RNMDT_relaxation = RNMDT_based_problem_generation(initial_parameters, generated_parameters)
[ fix.(RNMDT_relaxation[:x][:,s], x_k[s]) for s = 1 : initial_parameters.num_scen ]
print(RNMDT_relaxation)

# optimising the problem
optimize!(RNMDT_relaxation)

# plot the results
print("\n******************** RNMDT relaxation *******************\n\n")
print("objective function value : $(objective_value(RNMDT_relaxation))\n")
print("first stage variables : $(value.(RNMDT_relaxation[:x]))\n")
print("second stage variables : $(value.(RNMDT_relaxation[:y]))\n")
print("dual gap : $(MOI.get(RNMDT_relaxation, MOI.RelativeGap()))\n")
print("\n*********************************************************\n")


optimize!(try_node.dual_subproblems[2])
value.(try_node.dual_subproblems[2][:x])
value.(try_node.dual_subproblems[2][:y])
value.(try_node.dual_subproblems[2][:w_RNMDT])

@constraint(try_node.dual_subproblems[1],tighten_bounds[1][1]<=try_node.dual_subproblems[1][:y][1]<=tighten_bounds[1][2])
@constraint(try_node.dual_subproblems[2],tighten_bounds[2][1]<=try_node.dual_subproblems[2][:y][1]<=tighten_bounds[2][2])
#@constraint(try_node.dual_subproblems[2], try_node.dual_subproblems[2][:y][1]<=3.19726787)
@constraint(try_node.dual_subproblems[3], try_node.dual_subproblems[3][:x][2]<=3.19726787)
@constraint(try_node.dual_subproblems[4], try_node.dual_subproblems[4][:x][4]<=3.32056313)
result = proximal_bundle_method(try_node, initial_centre_of_gravity, 1e-9)
al_problem, RNMDT_subproblems , initial_parameters, generated_parameters

tighten_bounds[1][2]

!
initial_parameters.al_is_used
