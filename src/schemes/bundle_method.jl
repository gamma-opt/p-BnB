# the function contains the implementation of the Lagrangian decomposition
# with the bundle method multipliers update applied to the mixed integer based
# relaxation of the original problem

function auxiliary_check(vector_of_lambda_lagrangian, current_teration, iterations_to_be_considered, num_scen, tolerance)
    result = zeros(num_scen-1, 1)
    result_tolerance = tolerance .* ones(num_scen-1, 1)
    for i = 1 : iterations_to_be_considered
        result = result .+ norm.(vector_of_lambda_lagrangian[current_teration - i + 1 , :] .- vector_of_lambda_lagrangian[current_teration - i, :])
    end

    return sum(result .> result_tolerance)/ (num_scen-1)
end

function bundle_method(bnb_node::node, max_number_of_iterations::Int, center_of_gravity_inital_value::Array{Any}, parallelised::String)


    number_of_the_serious_steps = 0

    # auxiliary function for Lagrangian lagrangian_multipliers_representing_variableltipliers update
    f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1 : end]) : - lambda_lagrangian[dec_index-1])

    # lagrangian relaxation variables for the x and y non anticipativity conditions written in the column, for each iteration
    vector_of_lambda_lagrangian = Array{Any}(undef, max_number_of_iterations, bnb_node.initial_parameters.num_scen - 1)
    #[ vector_of_lambda_lagrangian[1, i] = lagrangian_multipliers_min .+ (lagrangian_multipliers_max - lagrangian_multipliers_min)  .* rand(1,
        #bnb_node.initial_parameters.num_int_var)
            #for i = 1 : bnb_node.initial_parameters.num_scen - 1 ]
    [ vector_of_lambda_lagrangian[1, i] = center_of_gravity_inital_value[i] for i = 1 : bnb_node.initial_parameters.num_scen - 1 ]

    # dual function at the lagragian multiplers' vector at correspondent iteration
    dual_objective_value_at_lagrangian = Array{Float64}(undef, max_number_of_iterations)

    # vector that contains decision variables written in a column (x and y in this case)
    # (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
    integer_decision_variables_values_for_each_scenario = Array{Float64}(undef, bnb_node.initial_parameters.num_int_var, bnb_node.initial_parameters.num_scen)

    # vector that contains decision variables written in a column (x and y in this case)
    # (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
    continuous_decision_variables_values_for_each_scenario = Array{Float64}(undef, bnb_node.initial_parameters.num_cont_var, bnb_node.initial_parameters.num_scen)


    # values at each iteration of the variable z uder minimization in the objective function of the cutting plane method
    relaxed_dual_objective_value = Array{Float64}(undef, 1, max_number_of_iterations)

    # the center of gravity at each iteration
    center_of_gravity = Array{Any}(undef, max_number_of_iterations, bnb_node.initial_parameters.num_scen - 1)

    # dual function at the center of gravity at correspondent iteration
    dual_function_value_at_the_center_of_gravity = Array{Float64}(undef, 1, max_number_of_iterations)

    # subgradient vector at each iteration
    subgradient_vector = Array{Any}(undef, max_number_of_iterations, bnb_node.initial_parameters.num_scen - 1)

    # upper bound for the original problem obtained by solving relaxed subproblems
    # and summing up the values of the objective functions
    ub_of_original_problem = Array{Float64}(undef,1,1)

    cutting_plane_subproblem = Model(with_optimizer(Gurobi.Optimizer, Threads = 1, Method = 4)) #, LogFile = loglink_par_bundle * "$(bnb_node.initial_parameters.num_scen)_scenarios_$(bnb_node.initial_parameters.num_cont_var)_cont_var_$(bnb_node.initial_parameters.num_int_var)_int_var_$(number_of_constraints)_constraints_$(seed)_seed_$(Dates.today())_bundle_LD+RNDMT_par_logfile.txt" ))
    @variables cutting_plane_subproblem begin
        z
        lagrangian_multipliers_representing_variable[ 1 : bnb_node.initial_parameters.num_int_var,
            1 : bnb_node.initial_parameters.num_scen - 1]
    end

    # the values for the parameters of the Bundle method
    m = 0.7
    d = 10

    iteration = 0 # strating counter

    # stopping criteria for the bundle method (percentage)
    eps_stop = 10 # stoping percentage
    number_of_iteration_for_checking = 3 # number of iterations
                                         # used to check the stopping criteria

    initial_time = time()
    #while (iteration < max_number_of_iterations) & ((iteration > number_of_iteration_for_checking + 1 ) ? ( norm((dual_objective_value_at_lagrangian[iteration - number_of_iteration_for_checking-1 : iteration-2] .- dual_objective_value_at_lagrangian[iteration - number_of_iteration_for_checking : iteration - 1]) ./ dual_objective_value_at_lagrangian[iteration - number_of_iteration_for_checking-1 : iteration-2] .* ones(number_of_iteration_for_checking) .*100 ) >= eps_stop) : true)
    while (iteration < max_number_of_iterations) & ((iteration > number_of_iteration_for_checking + 1 ) ? ( norm(dual_objective_value_at_lagrangian[iteration - number_of_iteration_for_checking : iteration-1] .- dual_objective_value_at_lagrangian[iteration - number_of_iteration_for_checking - 1 : iteration - 2]) >= eps_stop) : true)
        iteration += 1
        ub_of_original_problem[1] = 0

        if parallelised  == "parallelised"

                @suppress @sync Threads.@threads for s in 1 : bnb_node.initial_parameters.num_scen
                    #objective_update
                    @objective( bnb_node.dual_subproblems[s], Max,
                        bnb_node.initial_parameters.scen_prob[s] *
                        (
                        sum( bnb_node.dual_subproblems[s][:y][i] * bnb_node.generated_parameters.objective_Qs[s][i, j] * bnb_node.dual_subproblems[s][:y][j] for i = 1 : bnb_node.initial_parameters.num_cont_var, j = 1 : bnb_node.initial_parameters.num_cont_var)
                        + sum( bnb_node.dual_subproblems[s][:x][i] * bnb_node.generated_parameters.objective_c[i]  for i = 1:bnb_node.initial_parameters.num_int_var)
                        + sum( bnb_node.dual_subproblems[s][:y][j] * bnb_node.generated_parameters.objective_fs[s][j]  for j = 1:bnb_node.initial_parameters.num_cont_var)
                        )
                        +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s ) .* bnb_node.dual_subproblems[s][:x] )
                    )

                    #storing the optimal value of the obejective and variables of the p-LR with fixed multipliers
                    status = optimize!(bnb_node.dual_subproblems[s])
                    obj_value = objective_value(bnb_node.dual_subproblems[s])

                    ub_of_original_problem[1] = ub_of_original_problem[1] + obj_value
                    integer_decision_variables_values_for_each_scenario[ :, s ] = value.(bnb_node.dual_subproblems[s][:x])
                    continuous_decision_variables_values_for_each_scenario[ :, s ] = value.(bnb_node.dual_subproblems[s][:y])

            end

        else
            @suppress for s in 1 : bnb_node.initial_parameters.num_scen
                #objective_update
                @objective( bnb_node.dual_subproblems[s], Max,
                    bnb_node.initial_parameters.scen_prob[s] *
                    (
                    sum( bnb_node.dual_subproblems[s][:y][i] * bnb_node.generated_parameters.objective_Qs[s][i, j] * bnb_node.dual_subproblems[s][:y][j] for i = 1 : bnb_node.initial_parameters.num_cont_var, j = 1 : bnb_node.initial_parameters.num_cont_var)
                    + sum( bnb_node.dual_subproblems[s][:x][i] * bnb_node.generated_parameters.objective_c[i]  for i = 1:bnb_node.initial_parameters.num_int_var)
                    + sum(bnb_node.dual_subproblems[s][:y][j] * bnb_node.generated_parameters.objective_fs[s][j]  for j = 1:bnb_node.initial_parameters.num_cont_var)
                    )
                    +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s ) .* bnb_node.dual_subproblems[s][:x] )
                )

                #storing the optimal value of the obejective and variables of the p-LR with fixed multipliers
                status = optimize!(bnb_node.dual_subproblems[s])
                obj_value = objective_value(bnb_node.dual_subproblems[s])

                ub_of_original_problem[1] = ub_of_original_problem[1] + obj_value
                integer_decision_variables_values_for_each_scenario[ :, s ] = value.(bnb_node.dual_subproblems[s][:x])
                continuous_decision_variables_values_for_each_scenario[ :, s ] = value.(bnb_node.dual_subproblems[s][:y])

            end
        end


        dual_objective_value_at_lagrangian[iteration] = ub_of_original_problem[1]

        #calculating the subgradient
        [ subgradient_vector[iteration,  s - 1] =  integer_decision_variables_values_for_each_scenario[ :, 1] - integer_decision_variables_values_for_each_scenario[ :, s] for  s in 2: bnb_node.initial_parameters.num_scen ]

        if iteration == 1
            # if it is the first iteration simply set the centre of gravity to the initial values of the Lagrngian multipliers
            center_of_gravity[iteration, :] = vector_of_lambda_lagrangian[iteration, :]
            dual_function_value_at_the_center_of_gravity[iteration] = dual_objective_value_at_lagrangian[iteration]
    elseif  dual_function_value_at_the_center_of_gravity[iteration-1] - dual_objective_value_at_lagrangian[iteration] >= m * ( dual_function_value_at_the_center_of_gravity[iteration-1] -  ( relaxed_dual_objective_value[iteration-1] + d * sum( sum( (vector_of_lambda_lagrangian[iteration, s] .- center_of_gravity[iteration-1, s]) .^ 2 )  for s = 1 : bnb_node.initial_parameters.num_scen - 1 ) ) )
            center_of_gravity[iteration, :] = vector_of_lambda_lagrangian[iteration, :]
            dual_function_value_at_the_center_of_gravity[iteration] = dual_objective_value_at_lagrangian[iteration]
            number_of_the_serious_steps = number_of_the_serious_steps + 1
        else
            center_of_gravity[iteration, :] = center_of_gravity[iteration-1, :]
            dual_function_value_at_the_center_of_gravity[iteration] = dual_function_value_at_the_center_of_gravity[iteration-1]

        end

        @objective(cutting_plane_subproblem, Min, cutting_plane_subproblem[:z] + d * sum( sum( (cutting_plane_subproblem[:lagrangian_multipliers_representing_variable][:, s] .- center_of_gravity[iteration, s] ).^2 ) for  s in 1 : bnb_node.initial_parameters.num_scen - 1 ) )

        @constraint(cutting_plane_subproblem, cutting_plane_subproblem[:z] >= dual_function_value_at_the_center_of_gravity[iteration] + sum( sum( subgradient_vector[iteration, s] .* ( cutting_plane_subproblem[:lagrangian_multipliers_representing_variable][:, s] .- vector_of_lambda_lagrangian[iteration, s] ) ) for s = 1 : bnb_node.initial_parameters.num_scen - 1) )

        status = optimize!(cutting_plane_subproblem)

            if iteration < max_number_of_iterations
                [ vector_of_lambda_lagrangian[iteration+1, s] = value.(cutting_plane_subproblem[:lagrangian_multipliers_representing_variable][:, s]) for s = 1 : bnb_node.initial_parameters.num_scen - 1 ]
                relaxed_dual_objective_value[iteration] = value.(cutting_plane_subproblem[:z])
            end

    end

    final_time = time()-initial_time

    #return dual_objective_value_at_lagrangian[1:iteration], integer_decision_variables_values_for_each_scenario[:, :], continuous_decision_variables_values_for_each_scenario[:, :], vector_of_lambda_lagrangian[ 1 : iteration, :], number_of_the_serious_steps, center_of_gravity[iteration, :]
    return bundle_method_output(dual_objective_value_at_lagrangian[1:iteration], integer_decision_variables_values_for_each_scenario, continuous_decision_variables_values_for_each_scenario)
end
