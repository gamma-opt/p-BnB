# definining initial parameters 

# the number of scenarios
n_scenarios = 5
# the number of the first stage variables / per scenario
fs_var = 5
# the number of the second stage variabes / per scenario
ss_var = 5
# the number of the constraints / per scenario
const_num = 10
# RNMDT precision factor p
precision_factor = -1

# defining general intial parameters for the optimisation problem
initial_parameters = initialisation(n_scenarios, fs_var, ss_var,const_num, precision_factor)

# generating the structure containing the constraints and objective related parameters
generated_parameters = parameters_generation(initial_parameters) 

# defining the dual variables values 
w_0 = intial_centre_of_gravity_generation(initial_parameters.bm_parameters.initial_centre_of_gravity_min, initial_parameters.bm_parameters.initial_centre_of_gravity_max, initial_parameters.random_seed, (initial_parameters.al_is_used ? initial_parameters.num_scen : initial_parameters.num_scen - 1), initial_parameters.num_first_stage_var)

# formulating the agumented lagrangian problem
ag_lagranigan_problem  = RNMDT_based_augmented_lagrangian_relaxation_problem_generation(initial_parameters, generated_parameters, initial_parameters.RNMDT_precision_factor)

# defining the bnb_structure to call the functions
bnb_struct = bnb_model(initial_parameters)

# gethering the first node
current_node, current_node_id = node_selection(bnb_struct)

# formulating initial set V and x0 
V_0, x_0 = FW_PH_V_0_initialisation(current_node, w_0)

# formulating z_0 
z_0 = sum(initial_parameters.scen_prob[s] .* x_0[s]  for s = 1:initial_parameters.num_scen)

# folrmulating dual multipliers to eneter the SDM and satisfy dual feasibility condition
w_1 =  Array{Array{Float64}}(undef, initial_parameters.num_scen)
[ w_1[s] = w_0[s] .+ initial_parameters.al_penalty_parameter .* (x_0[s] .- z_0) for s = 1:initial_parameters.num_scen ]

# checking dual feasibility condition 
dual_feasibility_condition = zeros(length(w_1[1]))
for s = 1 : initial_parameters.num_scen
    dual_feasibility_condition = dual_feasibility_condition .+ initial_parameters.scen_prob[s] .* w_1[s]
end
println("dual feasibility condition: $(dual_feasibility_condition)")

##  Solving the augmented lagrangian problem with fixed values for dual variables and z

# defining the array to store the dual objective values of the augmented lagrangian scenario-subproblems solved with Gurobi
Gurobi_solved_ag_subproblems = []

# fixing the values of the variable z and dual multipliers in the objective of the augmented lagrangian 
for s = 1:initial_parameters.num_scen

    #defining the objective with the fixed values of the lagrangian multipliers
    @objective(ag_lagranigan_problem[s], Min,
    -
    ( sum(generated_parameters.objective_Qs[s][i, j] * ag_lagranigan_problem[s][:w_RNMDT][i, j]
        for i = 1 : initial_parameters.num_second_stage_var,
            j = 1 : initial_parameters.num_second_stage_var)
    + sum( ag_lagranigan_problem[s][:x][i] * generated_parameters.objective_c[i]  for i = 1:initial_parameters.num_first_stage_var)
    + sum( ag_lagranigan_problem[s][:y][j] * generated_parameters.objective_fs[s][j]  for j = 1:initial_parameters.num_second_stage_var)
    )

    + sum( w_1[s] .* (ag_lagranigan_problem[s][:x]) )
    + sum( initial_parameters.al_penalty_parameter/2 .* (ag_lagranigan_problem[s][:x] .- ag_lagranigan_problem[s][:al_z]) .* (ag_lagranigan_problem[s][:x] .- ag_lagranigan_problem[s][:al_z]) )

    + initial_parameters.μ * sum(ag_lagranigan_problem[s][:z][r] for r  = 1:initial_parameters.num_const )

    )
        
    # fixing the value of the variable al_z
    
    if sum(is_integer.((ag_lagranigan_problem[s][:al_z]))) == length((ag_lagranigan_problem[s][:al_z]))
        JuMP.unset_integer.(ag_lagranigan_problem[s][:al_z])
    end

    fix.(ag_lagranigan_problem[s][:al_z], z_0)

    optimize!(ag_lagranigan_problem[s])
    push!(Gurobi_solved_ag_subproblems, objective_value(ag_lagranigan_problem[s]))

end
Gurobi_solved_ag = sum( initial_parameters.scen_prob .* Gurobi_solved_ag_subproblems)
println("Gurobi solved instance: $Gurobi_solved_ag")


## solving the the agumented lagrangian with SDM 

# defining the array to store scenario-based dual objective values resulting from SDM 
SDM_solved_ag_subproblems = []

for s = 1:initial_parameters.num_scen
    SDM_output = SDM(s, current_node, V_0[s], x_0[s], V_0[s][1][2], V_0[s][1][3], V_0[s][1][4], w_1[s], z_0, 1000, 1e-8)
    push!(SDM_solved_ag_subproblems, SDM_output[end - 1 ])
end

SDM_solved_ag = sum(initial_parameters.scen_prob .* SDM_solved_ag_subproblems)

println("SDM solved isntance $SDM_solved_ag ")