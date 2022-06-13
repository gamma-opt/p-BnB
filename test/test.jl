include("initialization.jl")
primal_problem = MIP_generation(initial_parameters, generated_parameters)

optimize!(primal_problem)

# plot the results
print("\n********************* PRIMAL PROBLEM ********************\n\n")
print("objective function value : $(objective_value(primal_problem))\n")
print("first stage variables : $(value.(primal_problem[:x]))\n")
print("second stage variables : $(value.(primal_problem[:y]))\n")
print("dual gap : $(MOI.get(primal_problem, MOI.RelativeGap()))\n")
print("\n**********************************************************\n")

##
RNMDT_relaxation = RNMDT_based_problem_generation(initial_parameters, generated_parameters)
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
