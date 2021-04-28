function intial_centre_of_gravity_generation(centre_of_gravity_min, centre_of_gravity_max, random_seed, num_scen, num_first_stage_var)
    Random.seed!(random_seed)
    centre_of_gravity_inital_value = Array{Array{Float64}}(undef, num_scen)
    [centre_of_gravity_inital_value[i] = centre_of_gravity_min .+ (centre_of_gravity_max - centre_of_gravity_min) .* rand(num_first_stage_var)
        for i = 1 : num_scen]

    return centre_of_gravity_inital_value
end


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

function bundle_method(bnb_node::node, centre_of_gravity_inital_value)

    # simplifying the references
    input_parameters = bnb_node.initial_parameters.bm_parameters

    number_of_the_serious_steps = 0

    # auxiliary function for Lagrangian lagrangian_multipliers_representing_variableltipliers update
    f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1 : end]) : - lambda_lagrangian[dec_index-1])

    # lagrangian relaxation variables for the x and y non anticipativity conditions written in the column, for each iteration
     vector_of_lambda_lagrangian = Array{Any}(undef, input_parameters.max_number_of_iterations, bnb_node.initial_parameters.num_scen - 1)
    [ vector_of_lambda_lagrangian[1, i] = centre_of_gravity_inital_value[i] for i = 1 : bnb_node.initial_parameters.num_scen - 1 ]

    # dual function at the lagragian multiplers' vector at correspondent iteration
    dual_objective_value_at_lagrangian = Array{Float64}(undef, input_parameters.max_number_of_iterations)

    # vector that contains decision variables written in a column (x and y in this case)
    # (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
    first_stage_decision_variables_values_for_each_scenario = Array{Float64}(undef, bnb_node.initial_parameters.num_first_stage_var, bnb_node.initial_parameters.num_scen)

    # vector that contains decision variables written in a column (x and y in this case)
    # (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
    second_stage_decision_variables_values_for_each_scenario = Array{Float64}(undef, bnb_node.initial_parameters.num_second_stage_var, bnb_node.initial_parameters.num_scen)

    if  bnb_node.initial_parameters.RNMDT_is_used
        RNMDT_quadraticity_variables_w_for_each_scenario =  Array{Float64}(undef, bnb_node.initial_parameters.num_second_stage_var, bnb_node.initial_parameters.num_second_stage_var, bnb_node.initial_parameters.num_scen)
    end

    # values at each iteration of the variable z uder minimization in the objective function of the cutting plane method
    relaxed_dual_objective_value = Array{Float64}(undef, 1, input_parameters.max_number_of_iterations)

    # the center of gravity at each iteration
    center_of_gravity = Array{Any}(undef, input_parameters.max_number_of_iterations, bnb_node.initial_parameters.num_scen - 1)

    # dual function at the center of gravity at correspondent iteration
    dual_function_value_at_the_center_of_gravity = Array{Float64}(undef, 1, input_parameters.max_number_of_iterations)

    # subgradient vector at each iteration
    subgradient_vector = Array{Any}(undef, input_parameters.max_number_of_iterations, bnb_node.initial_parameters.num_scen - 1)

    # upper bound for the original problem obtained by solving relaxed subproblems
    # and summing up the values of the objective functions
    ub_of_original_problem = Array{Float64}(undef,1,1)

    cutting_plane_subproblem = Model(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" =>  bnb_node.initial_parameters.gurobi_parameters.MIPGap, "Method" => bnb_node.initial_parameters.gurobi_parameters.Method, "OutputFlag" => bnb_node.initial_parameters.gurobi_parameters.OutputFlag,  "Threads" => bnb_node.initial_parameters.gurobi_parameters.Threads)) #, LogFile = loglink_par_bundle * "$(bnb_node.initial_parameters.num_scen)_scenarios_$(bnb_node.initial_parameters.num_second_stage_var)_cont_var_$(bnb_node.initial_parameters.num_first_stage_var)_int_var_$(number_of_constraints)_constraints_$(seed)_seed_$(Dates.today())_bundle_LD+RNDMT_par_logfile.txt" ))
    @variables cutting_plane_subproblem begin
        z
        lagrangian_multipliers_representing_variable[ 1 : bnb_node.initial_parameters.num_first_stage_var,
            1 : bnb_node.initial_parameters.num_scen - 1]
    end

    iteration = 0 # strating counter


    initial_time = time()
    #while (iteration < input_parameters.max_number_of_iterations) & ((iteration > input_parameters.number_of_iteration_for_checking + 1 ) ? ( norm((dual_objective_value_at_lagrangian[iteration - input_parameters.number_of_iteration_for_checking-1 : iteration-2] .- dual_objective_value_at_lagrangian[iteration - input_parameters.number_of_iteration_for_checking : iteration - 1]) ./ dual_objective_value_at_lagrangian[iteration - input_parameters.number_of_iteration_for_checking-1 : iteration-2] .* ones(input_parameters.number_of_iteration_for_checking) .*100 ) >= input_parameters.eps_stop) : true)
    while (iteration < input_parameters.max_number_of_iterations) & ((iteration > input_parameters.number_of_iteration_for_checking + 1 ) ? ( norm(dual_objective_value_at_lagrangian[iteration - input_parameters.number_of_iteration_for_checking : iteration-1] .- dual_objective_value_at_lagrangian[iteration - input_parameters.number_of_iteration_for_checking - 1 : iteration - 2]) >= input_parameters.eps_stop) : true)
        iteration += 1
        ub_of_original_problem[1] = 0

        if input_parameters.parallelisation_is_used == true # if we use scenario-based parallelisation

                @suppress @sync Threads.@threads for s in 1 : bnb_node.initial_parameters.num_scen

                    # objective_update

                    if bnb_node.initial_parameters.RNMDT_is_used  # if we apply LR to the RNMDT-based relaxation

                        @objective( bnb_node.dual_subproblems[s], Max,
                            bnb_node.initial_parameters.scen_prob[s] *
                            (
                            sum( bnb_node.generated_parameters.objective_Qs[s][i, j] * bnb_node.dual_subproblems[s][:w_RNMDT][i,j]
                                for i = 1 : bnb_node.initial_parameters.num_second_stage_var,
                                    j = 1 : bnb_node.initial_parameters.num_second_stage_var)
                            + sum( bnb_node.dual_subproblems[s][:x][i] * bnb_node.generated_parameters.objective_c[i]  for i = 1:bnb_node.initial_parameters.num_first_stage_var)
                            + sum( bnb_node.dual_subproblems[s][:y][j] * bnb_node.generated_parameters.objective_fs[s][j]  for j = 1:bnb_node.initial_parameters.num_second_stage_var)
                            )
                            +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s ) .* bnb_node.dual_subproblems[s][:x] )
                        )

                    else # if we apply LR straight to the primal problem

                        @objective( bnb_node.dual_subproblems[s], Max,
                            bnb_node.initial_parameters.scen_prob[s] *
                            (
                            sum( bnb_node.dual_subproblems[s][:y][i] * bnb_node.generated_parameters.objective_Qs[s][i, j] * bnb_node.dual_subproblems[s][:y][j] for i = 1 : bnb_node.initial_parameters.num_second_stage_var, j = 1 : bnb_node.initial_parameters.num_second_stage_var)
                            + sum( bnb_node.dual_subproblems[s][:x][i] * bnb_node.generated_parameters.objective_c[i]  for i = 1:bnb_node.initial_parameters.num_first_stage_var)
                            + sum( bnb_node.dual_subproblems[s][:y][j] * bnb_node.generated_parameters.objective_fs[s][j]  for j = 1:bnb_node.initial_parameters.num_second_stage_var)
                            )
                            +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s ) .* bnb_node.dual_subproblems[s][:x] )
                        )

                    end

                    # storing the optimal value of the obejective and variables of the p-LR with fixed multipliers
                    status = optimize!(bnb_node.dual_subproblems[s])
                    obj_value = objective_value(bnb_node.dual_subproblems[s])

                    ub_of_original_problem[1] = ub_of_original_problem[1] + obj_value
                    first_stage_decision_variables_values_for_each_scenario[ :, s ] = value.(bnb_node.dual_subproblems[s][:x])
                    second_stage_decision_variables_values_for_each_scenario[ :, s ] = value.(bnb_node.dual_subproblems[s][:y])

                    if bnb_node.initial_parameters.RNMDT_is_used # if we apply LR to the RNMDT-based relaxation we store auxiliary variable w
                        RNMDT_quadraticity_variables_w_for_each_scenario[ :, :, s] = value.(bnb_node.dual_subproblems[s][:w_RNMDT])
                    end
            end # for loop

        else
            @suppress for s in 1 : bnb_node.initial_parameters.num_scen # if we don't use scenario-based parallelisation

                # objective_update

                if bnb_node.initial_parameters.RNMDT_is_used  # if we apply LR to the RNMDT-based relaxation

                    @objective( bnb_node.dual_subproblems[s], Max,
                        bnb_node.initial_parameters.scen_prob[s] *
                        (
                        sum( bnb_node.generated_parameters.objective_Qs[s][i, j] * bnb_node.dual_subproblems[s][:w_RNMDT][i,j]
                            for i = 1 : bnb_node.initial_parameters.num_second_stage_var,
                                j = 1 : bnb_node.initial_parameters.num_second_stage_var)
                        + sum( bnb_node.dual_subproblems[s][:x][i] * bnb_node.generated_parameters.objective_c[i]  for i = 1:bnb_node.initial_parameters.num_first_stage_var)
                        + sum( bnb_node.dual_subproblems[s][:y][j] * bnb_node.generated_parameters.objective_fs[s][j]  for j = 1:bnb_node.initial_parameters.num_second_stage_var)
                        )
                        +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s ) .* bnb_node.dual_subproblems[s][:x] )
                    )

                else # if we apply LR straight to the primal problem

                    @objective( bnb_node.dual_subproblems[s], Max,
                        bnb_node.initial_parameters.scen_prob[s] *
                        (
                        sum( bnb_node.dual_subproblems[s][:y][i] * bnb_node.generated_parameters.objective_Qs[s][i, j] * bnb_node.dual_subproblems[s][:y][j] for i = 1 : bnb_node.initial_parameters.num_second_stage_var, j = 1 : bnb_node.initial_parameters.num_second_stage_var)
                        + sum( bnb_node.dual_subproblems[s][:x][i] * bnb_node.generated_parameters.objective_c[i]  for i = 1:bnb_node.initial_parameters.num_first_stage_var)
                        + sum( bnb_node.dual_subproblems[s][:y][j] * bnb_node.generated_parameters.objective_fs[s][j]  for j = 1:bnb_node.initial_parameters.num_second_stage_var)
                        )
                        +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s ) .* bnb_node.dual_subproblems[s][:x] )
                    )

                end

                # storing the optimal value of the obejective and variables of the p-LR with fixed multipliers
                status = optimize!(bnb_node.dual_subproblems[s])
                obj_value = objective_value(bnb_node.dual_subproblems[s])

                ub_of_original_problem[1] = ub_of_original_problem[1] + obj_value
                first_stage_decision_variables_values_for_each_scenario[ :, s ] = value.(bnb_node.dual_subproblems[s][:x])
                second_stage_decision_variables_values_for_each_scenario[ :, s ] = value.(bnb_node.dual_subproblems[s][:y])

                if bnb_node.initial_parameters.RNMDT_is_used # if we apply LR to the RNMDT-based relaxation we store auxiliary variable w
                    RNMDT_quadraticity_variables_w_for_each_scenario[ :, :, s] = value.(bnb_node.dual_subproblems[s][:w_RNMDT])
                end

            end # for loop

        end # if (parallelisation based)


        dual_objective_value_at_lagrangian[iteration] = ub_of_original_problem[1]

        # calculating the subgradient
        [ subgradient_vector[iteration,  s - 1] =  first_stage_decision_variables_values_for_each_scenario[ :, 1] - first_stage_decision_variables_values_for_each_scenario[ :, s] for  s in 2 : bnb_node.initial_parameters.num_scen ]

        if iteration == 1
            # if it is the first iteration simply set the centre of gravity to the initial values of the Lagrngian multipliers
            center_of_gravity[iteration, :] = vector_of_lambda_lagrangian[iteration, :]
            dual_function_value_at_the_center_of_gravity[iteration] = dual_objective_value_at_lagrangian[iteration]
    elseif  dual_function_value_at_the_center_of_gravity[iteration-1] - dual_objective_value_at_lagrangian[iteration] >= input_parameters.m * ( dual_function_value_at_the_center_of_gravity[iteration-1] -  ( relaxed_dual_objective_value[iteration-1] + input_parameters.d * sum( sum( (vector_of_lambda_lagrangian[iteration, s] .- center_of_gravity[iteration-1, s]) .^ 2 )  for s = 1 : bnb_node.initial_parameters.num_scen - 1 ) ) )
            center_of_gravity[iteration, :] = vector_of_lambda_lagrangian[iteration, :]
            dual_function_value_at_the_center_of_gravity[iteration] = dual_objective_value_at_lagrangian[iteration]
            number_of_the_serious_steps = number_of_the_serious_steps + 1
        else
            center_of_gravity[iteration, :] = center_of_gravity[iteration-1, :]
            dual_function_value_at_the_center_of_gravity[iteration] = dual_function_value_at_the_center_of_gravity[iteration-1]

        end

        @objective(cutting_plane_subproblem, Min, cutting_plane_subproblem[:z] + input_parameters.d * sum( sum( (cutting_plane_subproblem[:lagrangian_multipliers_representing_variable][:, s] .- center_of_gravity[iteration, s] ).^2 ) for  s in 1 : bnb_node.initial_parameters.num_scen - 1 ) )

        @constraint(cutting_plane_subproblem, cutting_plane_subproblem[:z] >= dual_function_value_at_the_center_of_gravity[iteration] + sum( sum( subgradient_vector[iteration, s] .* ( cutting_plane_subproblem[:lagrangian_multipliers_representing_variable][:, s] .- vector_of_lambda_lagrangian[iteration, s] ) ) for s = 1 : bnb_node.initial_parameters.num_scen - 1) )

        @suppress status = optimize!(cutting_plane_subproblem)

            if iteration < input_parameters.max_number_of_iterations
                [ vector_of_lambda_lagrangian[iteration+1, s] = value.(cutting_plane_subproblem[:lagrangian_multipliers_representing_variable][:, s]) for s = 1 : bnb_node.initial_parameters.num_scen - 1 ]
                relaxed_dual_objective_value[iteration] = value.(cutting_plane_subproblem[:z])
            end

    end

    final_time = time()-initial_time

    print("number of the serious steps: $number_of_the_serious_steps\n")

    if bnb_node.initial_parameters.RNMDT_is_used # if we apply LR to the RNMDT-based relaxation we store auxiliary variable w
        return bm_output(dual_objective_value_at_lagrangian[1:iteration], [first_stage_decision_variables_values_for_each_scenario, second_stage_decision_variables_values_for_each_scenario, RNMDT_quadraticity_variables_w_for_each_scenario,  center_of_gravity[iteration]])
    else
        return bm_output(dual_objective_value_at_lagrangian[1:iteration], [first_stage_decision_variables_values_for_each_scenario, second_stage_decision_variables_values_for_each_scenario])
    end

end
