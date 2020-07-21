"""
    a (parent_node::Model, var_index::Int, inequality::String, value::Int)
returns a child node generated from a parent_node by copying it
and addinbg auxiliary constraint x[var_index] "<=" or ">=" value.

"""

function MIP_generation(intial_parameters::MIP_initial_parameters)

    # generating the parameters
    generated_parameters = parameters_generation(initial_parameters)

    # depending on whether the variables are fixed or not we use Ipopt and Gurobi as a solver respectively
    original_problem = Model( (initial_parameters.is_int_fixed == true) ? ( with_optimizer(Ipopt.Optimizer, print_level=0) ) : ( with_optimizer(Gurobi.Optimizer, NonConvex = 2, MIPGap =  0, Method = 4, OutputFlag=0, TimeLimit = initial_parameters.solver_time_limit, Threads = 1) ) )

    # integer decision variables
    if (initial_parameters.is_int_fixed == true)
        @variable(original_problem, x[ 1 : initial_parameters.num_int_var, 1 : initial_parameters.num_scen ] )
    else
        @variable(original_problem, x[ 1 : initial_parameters.num_int_var, 1 : initial_parameters.num_scen ], Int )
    end

    # continuous decision variables
    @variable(original_problem, y[ 1 : initial_parameters.num_cont_var, 1 : initial_parameters.num_scen ])

    # quadratic objective
    @objective(original_problem, Max,
        sum( initial_parameters.scen_prob[s] *
            (
                sum( y[i, s] * generated_parameters.objective_Qs[s][i, j] * y[j, s] for i = 1 : initial_parameters.num_cont_var, j = 1 : initial_parameters.num_cont_var)
                + sum( x[i, s] * generated_parameters.objective_c[i] for i = 1 : initial_parameters.num_int_var)
                + sum( y[i, s] * generated_parameters.objective_fs[s][i] for i = 1 : initial_parameters.num_cont_var)
            )
        for s in 1:initial_parameters.num_scen)
        )

    # quadratic constraints
    @constraint(original_problem, [s = 1 : initial_parameters.num_scen, i = 1 : initial_parameters.num_const ],
        sum( y[j, s] * generated_parameters.constraint_Qs[s][i][j,k] * y[k, s] for j = 1 : initial_parameters.num_cont_var, k = 1: initial_parameters.num_cont_var)
        + sum( x[j, s] * generated_parameters.constraint_fs[s][i][1, j] for j = 1 : initial_parameters.num_int_var)
        + sum( y[j, s] * generated_parameters.constraint_fs[s][i][2, j] for j = 1:  initial_parameters.num_cont_var)
        + generated_parameters.constraint_fs[s][i][3, 1] <= 0 )

    # box constraints for integer variables
    @constraint(original_problem, [ s = 1 : initial_parameters.num_scen ],
        generated_parameters.x_boundaries[:, 1] .<= x[:, s] .<= generated_parameters.x_boundaries[:, 2])

    # box constraints for continuous variables
    @constraint(original_problem, [s = 1 : initial_parameters.num_scen],
        generated_parameters.y_boundaries[:, 1] .<= y[:, s] .<= generated_parameters.y_boundaries[:, 2])

    # non-anticipativity conditions
    @constraint( original_problem, [s in 2 : initial_parameters.num_scen, i = 1 : initial_parameters.num_int_var],
        x[i, s] - x[i, 1] == 0 )

    return original_problem

end



# auxiliary function for Lagrangian multipliers update
f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1:end]) : - lambda_lagrangian[dec_index-1])

function MIP_lagrangian_relaxation_generation(initial_parameters::MIP_initial_parameters)

    # generating the parameters
    generated_parameters = parameters_generation(initial_parameters)


    # generate
    vector_of_lambda_lagrangian = Array{Any}(undef, initial_parameters.num_scen - 1)
    [ vector_of_lambda_lagrangian[i] = zeros(1,initial_parameters.num_int_var)
            for i = 1 : initial_parameters.num_scen - 1 ]

    # creating the array of subproblems
    vector_of_subproblems = Array{Any}(undef, initial_parameters.num_scen)

    # formulating the subproblems based on the scenarios
    for s = 1 : initial_parameters.num_scen
        vector_of_subproblems[s] = Model(with_optimizer(Gurobi.Optimizer, Method = 4, OutputFlag=0, MIPGap =  0, Threads = 1))

        # integer decision variables
        @variable( vector_of_subproblems[s], x[1 : initial_parameters.num_int_var], Int )

        # continuous decision variables
        @variable( vector_of_subproblems[s], y[1 : initial_parameters.num_cont_var] )


        # quadratic objective
        @objective( vector_of_subproblems[s], Max,
            initial_parameters.scen_prob[s] *
            (
            sum( y[i] * generated_parameters.objective_Qs[s][i, j] * y[j] for i = 1 : initial_parameters.num_cont_var, j = 1 : initial_parameters.num_cont_var)
            + sum( x[i] * generated_parameters.objective_c[i]  for i = 1:initial_parameters.num_int_var)
            + sum( y[j] * generated_parameters.objective_fs[s][j]  for j = 1:initial_parameters.num_cont_var)
            )

            +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[:], s ) .* x )

        )

        #quadratic constraints
        @constraint( vector_of_subproblems[s], [ r = 1 : initial_parameters.num_const ],
            sum( y[i] * generated_parameters.constraint_Qs[s][r][i, j] * y[j] for i = 1 : initial_parameters.num_cont_var, j = 1 : initial_parameters.num_cont_var)
            + sum( x[i] * generated_parameters.constraint_fs[s][r][1, i] for i = 1:initial_parameters.num_int_var)
            + sum( y[j] * generated_parameters.constraint_fs[s][r][2, j] for j = 1:initial_parameters.num_cont_var )
            + generated_parameters.constraint_fs[s][r][3, 1] <= 0
        )

        # box constraints for integer variables
        @constraint( vector_of_subproblems[s],
            generated_parameters.x_boundaries[:, 1] .<= x .<= generated_parameters.x_boundaries[:, 2])

        # box constraints for continuous variables
        @constraint( vector_of_subproblems[s],
            generated_parameters.y_boundaries[:, 1] .<= y .<= generated_parameters.y_boundaries[:, 2])

    end

    return vector_of_subproblems

end
