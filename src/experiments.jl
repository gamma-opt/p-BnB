#src_link  =  "/scratch/work/belyakn1/BnB_p_lagrangian/src/"
src_link  =  "/Users/nikitabelyak/Dropbox (Aalto)/branch-and-bound-caroe-and-schultz/src/"
#src_link =  "/Users/Fabricio/Documents/GitHub/BnB_Caroe-and-Schultz/src/"

#cd(chop(src_link, tail = 4))
using Pkg
Pkg.activate(".")
#Pkg.update()

Pkg.instantiate()

include(src_link*"initialization.jl")

# set unique envinronment for Gurobi
const GRB_ENV = Gurobi.Env()

# check whether the folder for the "today" experiments exists
# and if not create one 
if !isdir(chop(src_link, tail = 4) * "experiments_" * string(Dates.today()))
    mkdir(chop(src_link, tail = 4) * "experiments_" * string(Dates.today()))
end

output_link = chop(src_link, tail = 4) * "experiments_" * string(Dates.today()) * "/"

## Defining the initial parameters for the experiments

# the number of scenarios
n_scenarios = [20]
# the number of the first stage variables / per scenario
fs_var = [5]
# the number of the second stage variabes / per scenario
ss_var = [10]
# the number of the constraints / per scenario
const_num = [10]
# Methods to be used [full scale, RNMDT, BnB]: 1 = use, 0 = don't use
experiments_methods = [0, 0, 1]
# The maximum time limit for the set of instances
g_time_limit = 3600
# precision factor values to be considered for RNMDT
p_value = [-1]

experiments_function(n_scenarios, fs_var, ss_var, const_num, experiments_methods, g_time_limit, p_value, output_link)

experiments_function(2, 10, 10, 10, [0,0,1], g_time_limit, [-1], output_link)


##----------------------------HERE THE DEBUGGING CODE IS LOCATED-------------------------------------------------

## Generating and solving full-scale problem 

# defining general intial parameters for the optimisation problem
initial_parameters = initialisation(10,10,10,10,-2)

# generating the structure containing the constraints and objective related parameters
generated_parameters = parameters_generation(initial_parameters)


# soving the full scale problem using global solver
primal_problem = MIP_generation(initial_parameters, generated_parameters)

#print(primal_problem)
optimize!(primal_problem)
obj_2 = objective_value(primal_problem)

value.(primal_problem[:x])
value.(primal_problem[:z])
value.(primal_problem[:y])

# soving the full scale problem using global solver
RNMDT_problem = RNMDT_based_problem_generation(initial_parameters, generated_parameters)

#print(RNMDT_problem)
optimize!(RNMDT_problem)
objective_value(RNMDT_problem)
RNMDT_x = value.(RNMDT_problem[:x])
value.(RNMDT_problem[:z])
value.(RNMDT_problem[:y])

for s = 1:initial_parameters.num_scen
    for i = 1:initial_parameters.num_first_stage_var
        fix.(primal_problem[:x][i,s], RNMDT_x[i,s]) 
    end
end


model, bub = pooling_problem_one_layer_pools(chop(src_link, tail = 4) * "pooling_problem_data", [1.0])
print(model)

Variables, Objective_c, Affine_constraints, Quadratic_constraints, Box_constraints, bub = polling_problem_parameters_extraction(model, bub)
@show Quadratic_constraints[1].quadratic_coefficients_matrix



ip, gp, var, bub = polling_problem_parameters_generation(model, 4, 5, 5, -1, bub)
stoc_ppm = stochastic_pooling_problem_parameters(10.0, 10.0, 0.0, 1.0, 1.0, 1.0 )

# soving the full scale problem using global solver
primal_problem = pooling_MIP_generation(ip, gp, var, bub, stoc_ppm)
print(primal_problem)
optimize!(primal_problem)

@show objective_value(primal_problem)
print("Primal problem \n")
print("X: \n")
@show value.(primal_problem[:x])
print("Y: \n")
@show value.(primal_problem[:y])


RNMDT_problem = RNMDT_based_pooling_problem_generation(ip, gp, var, bub, stoc_ppm)
print(RNMDT_problem)
optimize!(RNMDT_problem)
@show objective_value(RNMDT_problem)
@show value.(RNMDT_problem[:y])
@show value.(RNMDT_problem[:x])
y = []

for p = -1:-1:-100
    @suppress global ip, gp, var, bub = polling_problem_parameters_generation(model, 2, 10, 10, p, bub)
    @suppress stoc_ppm = stochastic_pooling_problem_parameters(10.0, 10.0, 0.0, 5.0, 1.0, 3.0 )
    @suppress RNMDT_problem = RNMDT_based_pooling_problem_generation(ip, gp, var, bub, stoc_ppm)
    @suppress optimize!(RNMDT_problem)
    @show objective_value(RNMDT_problem)
    #@show value.(RNMDT_problem[:x])
    y = value.(RNMDT_problem[:y])
    #print(RNMDT_problem)
end

RNMDT_problem = RNMDT_based_pooling_problem_generation(ip, gp, var, bub, stoc_ppm)
print(RNMDT_problem)
optimize!(RNMDT_problem)
objective_value(RNMDT_problem)
print("RNMDT problem \n")
print("X: \n")
value.(RNMDT_problem[:x])
print("Y: \n")
value.(RNMDT_problem[:y])

al_subproblems = RNMDT_based_augmented_lagrangian_relaxation_pooling_problem_generation(ip, gp, var, bub, stoc_ppm)
print((al_subproblems[1]))
optimize!(al_subproblems[1])
al_1 = objective_value(al_subproblems[1])