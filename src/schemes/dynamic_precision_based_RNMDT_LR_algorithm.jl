
function dynamic_precision_RNMDT_algorithm(bnb_node::node)

    # creating a copy of the bnb_node to work with it and not alter the original
    dp_RNMDT_node = node_copy(bnb_node)

    # simplifying the references
    initial_parameters = dp_RNMDT_node.initial_parameters
    dp_RNMDT_parameters = initial_parameters.dp_RNMDT_parameters
    generated_parameters = dp_RNMDT_node.generated_parameters

    # intialisation of missing parameters
    precision_p = -1 .* ones(initial_parameters.num_second_stage_var, initial_parameters.num_scen)
    iteration = 0 # iteration counter
    UB_optimal = Inf # initial value for the upper bound to enter the loop
    LB_optimal = -Inf # initial value for the lower bound to enter the loop

    # generating the value of the N1 parameter based on the given  percentage of the total number of the second stage decision variables in the primal problem
    N1 = Int(round( dp_RNMDT_parameters.N1_percentage/100*(initial_parameters.num_second_stage_var)*initial_parameters.num_scen, digits = 0 ) )

    # auxiliary variable for ranking the precion factor values
    f_rank = zeros(initial_parameters.num_second_stage_var, initial_parameters.num_scen)

    # generating the inital values for the center of gravity
    center_of_gravity_min = 0
    center_of_gravity_max = 0

    # generating centre of gravity intial values
    centre_of_gravity_inital_value = intial_centre_of_gravity_generation(initial_parameters.bm_parameters.initial_centre_of_gravity_min, initial_parameters.bm_parameters.initial_centre_of_gravity_max, initial_parameters.random_seed, initial_parameters.num_scen, initial_parameters.num_first_stage_var)

    bounds_results = zeros(3, dp_RNMDT_parameters.max_number_of_iterations)
    variables_values_results = Array{Array{Float64}}(undef, 2, dp_RNMDT_parameters.max_number_of_iterations)

    init_time  = time()
    while ( iteration == 0 ? true : ( (UB_optimal - LB_optimal)/LB_optimal > dp_RNMDT_parameters.tolerance ) ) & ( time() - init_time < dp_RNMDT_parameters.time_limit ) & (iteration < dp_RNMDT_parameters.max_number_of_iterations)

        iteration  = iteration + 1


        dp_RNMDT_node.dual_subproblems = RNMDT_based_lagrangian_relaxation_problem_generation(initial_parameters, generated_parameters, precision_p)

        @time bm_output = bundle_method(dp_RNMDT_node, centre_of_gravity_inital_value)
        UB_new = bm_output.dual_objective_value[end]
        xr_new = bm_output.variables_values[1]
        yr_new = bm_output.variables_values[2]
        if length(bm_output.variables_values) == 2
            error("the field 'RNMDT_is_used' is fixed to 'false' in initial parameters")
        else
            w_RNMDT_r_new = bm_output.variables_values[3]
            center_of_gravity_inital_value = bm_output.variables_values[4]
        end

        if (UB_new < UB_optimal)
            UB_optimal = UB_new
            xr_UB_optimal = xr_new
            yr_UB_optimal = yr_new
        end

        #xr = repeat( sum(x_values_LD_RNMDT, dims = 2)./initial_parameters.num_scen, 1, initial_parameters.num_scen )
        xr_new = repeat( Int.(round.( sum(xr_new, dims = 2)./initial_parameters.num_scen, digits = 0)), 1, initial_parameters.num_scen)

        variables_values_results[1,iteration] = xr_new
        variables_values_results[2,iteration] = yr_new

        # setting strating values for the continuous variables
        #set_start_value.(original_problem[:y], yr_new)

        #fixing the values for the integer variables
        fix.(dp_RNMDT_node.primal_problem[:x], xr_new)

        @time @suppress optimize!(dp_RNMDT_node.primal_problem)
        LB_new  = objective_value(dp_RNMDT_node.primal_problem)


        if (LB_new > LB_optimal)
            LB_optimal = LB_new
            xr_LB_optimal = value.(dp_RNMDT_node.primal_problem[:x])
            yr_LB_optimal = value.(dp_RNMDT_node.primal_problem[:y])
        end


        # if iteration is not a multiple of N2
        if mod(iteration+1, dp_RNMDT_parameters.N2) != 0
            [ f_rank[j, s] = sum( abs(generated_parameters.constraint_Qs[s][r][i, j] * (w_RNMDT_r_new[i, j, s] - yr_new[i, s] * yr_new[j, s]))  for i = 1 : initial_parameters.num_second_stage_var, r = 1 : initial_parameters.num_const) for j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen]
            # f_rank = reshape(f_rank, 1, initial_parameters.num_second_stage_var + initial_parameters.num_scen)
            # sorted_N1_indeces = sortperm(f_rank[1 : end], rev = true)[1 : min(N1, initial_parameters.num_second_stage_var)]
            for i  = 1 : min(N1, initial_parameters.num_second_stage_var * initial_parameters.num_scen)
                precision_p[findall(x->x==maximum(f_rank),f_rank)[1]] = precision_p[findall(x->x==maximum(f_rank),f_rank)[1]] - 1
                f_rank[findall(x->x==maximum(f_rank),f_rank)[1]] = - Inf
            end

        else
            # find indexes with the biggest p
            indixes_max_p = findall(x -> x == maximum(precision_p), precision_p)
            precision_p[indixes_max_p] = precision_p[indixes_max_p] .- 1

        end




    bounds_results[:, iteration] = [UB_optimal, LB_optimal, (UB_optimal - LB_optimal)/LB_optimal]

    end

    return dp_based_RNMDT_output(bounds_results[1,1:iteration], [variables_values_results[1, iteration], variables_values_results[2, iteration]])

end
