"""
    a (parent_node::Model, var_index::Int, inequality::String, value::Int)
returns a child node generated from a parent_node by copying it
and addinbg auxiliary constraint x[var_index] "<=" or ">=" value.

"""

function MIP_generation(intial_parameters::MIP_inital_parameters)

    # generating the parameters
    generated_parameters = parameters_generation(initial_parameters)

    # depending on whether the variables are fixed or not we use Ipopt and Gurobi as a solver respectively
    original_problem = Model( (initial_parameters.is_int_fixed == true) ? ( with_optimizer(Ipopt.Optimizer, print_level=0) ) : ( with_optimizer(Gurobi.Optimizer, NonConvex = 2, MIPGap =  0, Method = 4, OutputFlag=0, TimeLimit = inital_parameters.solver_time_limit, Threads = 1) ) )

    # integer decision variables
    if (initial_parameters.is_int_fixed == true)
        @variable(original_problem, x[ 1 : inital_parameters.num_int_var, 1 : inital_parameters.num_scen ] )
    else
        @variable(original_problem, x[ 1 : inital_parameters.num_int_var, 1 : inital_parameters.num_scen ], Int )
    end

    # continuous decision variables
    @variable(original_problem, y[ 1 : inital_parameters.num_cont_var, 1 : inital_parameters.num_scen ])

    # quadratic objective
    @objective(original_problem, Max,
        sum( inital_parameters.scen_prob[s] *
            (
                sum( y[i, s] * generated_parameters.objective_Qs[s][i, j] * y[j, s] for i = 1 : inital_parameters.num_cont_var, j = 1 : inital_parameters.num_cont_var)
                + sum( x[i, s] * generated_parameters.objective_c[i] for i = 1 : inital_parameters.num_int_var)
                + sum( y[i, s] * generated_parameters.objective_fs[s][i] for i = 1 : inital_parameters.num_cont_var)
            )
        for s in 1:inital_parameters.num_scen)
        )

    # quadratic constraints
    @constraint(original_problem, [s = 1 : inital_parameters.num_scen, i = 1 : inital_parameters.num_const ],
        sum( y[j, s] * generated_parameters.constraint_Qs[s][i][j,k] * y[k, s] for j = 1 : inital_parameters.num_cont_var, k = 1: inital_parameters.num_cont_var)
        + sum( x[j, s] * generated_parameters.constraint_fs[s][i][1, j] for j = 1 : inital_parameters.num_int_var)
        + sum( y[j, s] * generated_parameters.constraint_fs[s][i][2, j] for j = 1:  inital_parameters.num_cont_var)
        + generated_parameters.constraint_fs[s][i][3, 1] <= 0 )

    # box constraints for integer variables
    @constraint(original_problem, [ s = 1 : inital_parameters.num_scen ],
        generated_parameters.x_boundaries[:, 1] .<= x[:, s] .<= generated_parameters.x_boundaries[:, 2])

    # box constraints for continuous variables
    @constraint(original_problem, [s = 1 : inital_parameters.num_scen],
        generated_parameters.y_boundaries[:, 1] .<= y[:, s] .<= generated_parameters.y_boundaries[:, 2])

    # non-anticipativity conditions
    @constraint( original_problem, [s in 2 : inital_parameters.num_scen, i = 1 : inital_parameters.num_int_var],
        x[i, s] - x[i, 1] == 0 )

    return original_problem

end



# auxiliary function for Lagrangian multipliers update
f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1:end]) : - lambda_lagrangian[dec_index-1])

function MIP_lagrangian_relaxation_generation(initial_parameters::MIP_inital_parameters)

    # generating the parameters
    generated_parameters = parameters_generation(initial_parameters)

    # randomized lagrangian relaxation variables for the x and y non anticipativity conditions written in the column
    Random.seed!(seed)
    vector_of_lambda_lagrangian = Array{Any}(undef, number_of_scenarios - 1)
    [ vector_of_lambda_lagrangian[i] = 0 .+ 0.0 .* rand(1, number_of_integer_decision_variables)
            for i = 1 : number_of_scenarios - 1 ]

    # creating the array of subproblems
    subproblem = Array{Any}(undef, 1, number_of_scenarios)

    # formulating the subproblems based on the scenarios
    for s = 1 : number_of_scenarios

        subproblem[s] = Model(with_optimizer(Gurobi.Optimizer, Method = 4, OutputFlag=0, MIPGap =  0, Threads = 1))

        # integer decision variables
        @variable( subproblem[s], x[1 : number_of_integer_decision_variables], Int )

        # continuous decision variables
        @variable( subproblem[s], y[1 : number_of_continuous_decision_variables] )

        # RNMDT variables
        @variables subproblem[s] begin
            y_heat_RNMDT[ 1 : number_of_continuous_decision_variables, 1 : number_of_continuous_decision_variables, minimum(precision_p) : -1 ] # x_heat
            delta_y_RNMDT[ 1 : number_of_continuous_decision_variables ] # delta_y_heat
            z_RNMDT[ 1 : number_of_continuous_decision_variables, minimum(precision_p) : -1 ], Bin
            w_RNMDT[ 1 : number_of_continuous_decision_variables, 1 : number_of_continuous_decision_variables ]
            delta_w_RNMDT[1 : number_of_continuous_decision_variables, 1 : number_of_continuous_decision_variables ]

        end

        # quadratic objective
        @objective( subproblem[s], Max,
            round( (1/number_of_scenarios), digits = 3) *

            ( sum(objective_Qs[s][i, j] * w_RNMDT[i, j]
                for i = 1 : number_of_continuous_decision_variables,
                    j = 1 : number_of_continuous_decision_variables)
            + sum( x[i] * objective_c[i]  for i = 1:number_of_integer_decision_variables)
            + sum( y[j] * objective_fs[s][j]  for j = 1:number_of_continuous_decision_variables)
            )

            +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[:], s ) .* x )

        )

        #quadratic constraint
        @constraint( subproblem[s], [ r = 1 : number_of_constraints ],
            sum(constraint_Qs[s][r][i, j] * w_RNMDT[i, j]
                for i = 1 : number_of_continuous_decision_variables,
                j = 1 : number_of_continuous_decision_variables)
            + sum( x[i] * constraint_fs[s][r][1, i] for i = 1:number_of_integer_decision_variables)
            + sum( y[j] * constraint_fs[s][r][2, j] for j = 1:number_of_continuous_decision_variables )
            + constraint_fs[s][r][3, 1] <= 0
        ) # 27

## auxiliary RNMDT constraints

        @constraint( subproblem[s], [ j  = 1 : number_of_continuous_decision_variables],
            y[j] == (y_boundaries[j, 2] - y_boundaries[j, 1])  *
            ( sum( 2.0^l * z_RNMDT[j, l] for l = precision_p[j, s] : -1) + delta_y_RNMDT[j] ) ) # 28

        @constraint( subproblem[s], [ i  = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables ],
            w_RNMDT[i, j] == (y_boundaries[j, 2] - y_boundaries[j, 1]) *
            ( sum( 2.0^l * y_heat_RNMDT[i, j, l] for l = precision_p[j, s] : -1) + delta_w_RNMDT[i, j] ) ) # 29

        @constraint( subproblem[s], [ j  = 1 : number_of_continuous_decision_variables ],
            0 <= delta_y_RNMDT[j] <= 2.0^precision_p[j, s] ) # 30

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables],
            2.0^precision_p[i, s] * ( y[i] - y_boundaries[i, 2] ) + y_boundaries[i, 2] * delta_y_RNMDT[j] <= delta_w_RNMDT[i, j]) # 31

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables],
            delta_w_RNMDT[i, j] <= 2.0^precision_p[i, s] * ( y[i] - y_boundaries[i, 1] ) + y_boundaries[i, 1] * delta_y_RNMDT[j]) # 31

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables],
            y_boundaries[i, 1] * delta_y_RNMDT[j] <= delta_w_RNMDT[i, j]) #32

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables],
            delta_w_RNMDT[i, j] <= y_boundaries[i, 2] * delta_y_RNMDT[j]) # 32

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, l = precision_p[j, s] : -1],
            y_boundaries[i, 1] * z_RNMDT[j, l]  <= y_heat_RNMDT[i, j, l]) # 33

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, l = precision_p[j, s] : -1],
            y_heat_RNMDT[i, j, l] <= y_boundaries[i, 2] * z_RNMDT[j, l]) # 33

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, l = precision_p[j, s] : -1],
            y_boundaries[i, 1] * (1 - z_RNMDT[j, l]) <= y[i] - y_heat_RNMDT[i, j, l] ) # 34

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, l = precision_p[j, s] : -1],
            y[i] - y_heat_RNMDT[i, j, l] <= y_boundaries[i, 2] * (1 - z_RNMDT[j, l]) ) # 34
##
        # box constraints for integer variables
        @constraint( subproblem[s],
            x_boundaries[:, 1] .<= x .<= x_boundaries[:, 2])

        # box constraints for continuous variables
        @constraint( subproblem[s],
            y_boundaries[:, 1] .<= y .<= y_boundaries[:, 2])

    end

    return subproblem

end
