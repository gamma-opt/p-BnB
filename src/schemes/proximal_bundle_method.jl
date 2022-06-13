# the function calculates the values for the objective function
# of the dual subproblems using predefined values of the Lagrangian multipliers
function dual_subproblems_values(bnb_node::node, lm_values)

    # vector that contains decision variables written in a column (x and y in this case)
    # (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
    fs_dv_values = Array{Float64}(undef, bnb_node.initial_parameters.num_first_stage_var, bnb_node.initial_parameters.num_scen)

    # vector that contains decision variables written in a column (x and y in this case)
    # (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
    ss_dv_values = Array{Float64}(undef, bnb_node.initial_parameters.num_second_stage_var, bnb_node.initial_parameters.num_scen)

    # if we apply LR to the RNMDT-based relaxation we store auxiliary variable w values for each subproblem
    if  bnb_node.initial_parameters.RNMDT_is_used
        RNMDT_w_values =  Array{Float64}(undef, bnb_node.initial_parameters.num_second_stage_var, bnb_node.initial_parameters.num_second_stage_var, bnb_node.initial_parameters.num_scen)
    end

    # vector that contains the values for the auxiliary variales ensuring full recourse
    z_values = Array{Float64}(undef, bnb_node.initial_parameters.num_const, bnb_node.initial_parameters.num_scen)

    # defining the vector that will contain the values of the dual objective functions for each subproblem
    dual_subrpoblems_objective_values_vector = Array{Float64}(undef,bnb_node.initial_parameters.num_scen)

    if bnb_node.initial_parameters.bm_parameters.parallelisation_is_used == true # if we use scenario-based parallelisation

            @suppress @sync Threads.@threads for s in 1 : bnb_node.initial_parameters.num_scen

                # objective_update

                if bnb_node.initial_parameters.RNMDT_is_used  # if we apply LR to the RNMDT-based relaxation

                    if bnb_node.initial_parameters.al_is_used
                        @objective( bnb_node.dual_subproblems[s], Min,
                        - bnb_node.initial_parameters.scen_prob[s] *
                            (
                            sum( bnb_node.generated_parameters.objective_Qs[s][i, j] * bnb_node.dual_subproblems[s][:w_RNMDT][i,j]
                                for i = 1 : bnb_node.initial_parameters.num_second_stage_var,
                                    j = 1 : bnb_node.initial_parameters.num_second_stage_var)
                            + sum( bnb_node.dual_subproblems[s][:x][i] * bnb_node.generated_parameters.objective_c[i]  for i = 1:bnb_node.initial_parameters.num_first_stage_var)
                            + sum( bnb_node.dual_subproblems[s][:y][j] * bnb_node.generated_parameters.objective_fs[s][j]  for j = 1:bnb_node.initial_parameters.num_second_stage_var)
                            - bnb_node.initial_parameters.μ * sum(bnb_node.dual_subproblems[s][:z][r] for r  = 1:bnb_node.initial_parameters.num_const )
                            )
                            +  sum( lm_values[s] .* bnb_node.dual_subproblems[s][:x] )
                        )
                    else
                        @objective( bnb_node.dual_subproblems[s], Min,
                        - bnb_node.initial_parameters.scen_prob[s] *
                            (
                            sum( bnb_node.generated_parameters.objective_Qs[s][i, j] * bnb_node.dual_subproblems[s][:w_RNMDT][i,j]
                                for i = 1 : bnb_node.initial_parameters.num_second_stage_var,
                                    j = 1 : bnb_node.initial_parameters.num_second_stage_var)
                            + sum( bnb_node.dual_subproblems[s][:x][i] * bnb_node.generated_parameters.objective_c[i]  for i = 1:bnb_node.initial_parameters.num_first_stage_var)
                            + sum( bnb_node.dual_subproblems[s][:y][j] * bnb_node.generated_parameters.objective_fs[s][j]  for j = 1:bnb_node.initial_parameters.num_second_stage_var)
                            - bnb_node.initial_parameters.μ * sum(bnb_node.dual_subproblems[s][:z][r] for r  = 1:bnb_node.initial_parameters.num_const )
                            )
                            +  sum( f_lambda_lagrangian( lm_values[:], s ) .* bnb_node.dual_subproblems[s][:x] )
                        )
                    end

                else # if we apply LR straight to the primal problem

                    @objective( bnb_node.dual_subproblems[s], Min,
                      - bnb_node.initial_parameters.scen_prob[s] *
                        (
                        sum( bnb_node.dual_subproblems[s][:y][i] * bnb_node.generated_parameters.objective_Qs[s][i, j] * bnb_node.dual_subproblems[s][:y][j] for i = 1 : bnb_node.initial_parameters.num_second_stage_var, j = 1 : bnb_node.initial_parameters.num_second_stage_var)
                        + sum( bnb_node.dual_subproblems[s][:x][i] * bnb_node.generated_parameters.objective_c[i]  for i = 1:bnb_node.initial_parameters.num_first_stage_var)
                        + sum( bnb_node.dual_subproblems[s][:y][j] * bnb_node.generated_parameters.objective_fs[s][j]  for j = 1:bnb_node.initial_parameters.num_second_stage_var)
                        - bnb_node.initial_parameters.μ * sum(bnb_node.dual_subproblems[s][:z][r] for r  = 1:bnb_node.initial_parameters.num_const )
                        )
                        +  sum( f_lambda_lagrangian( lm_values[:], s ) .* bnb_node.dual_subproblems[s][:x] )
                    )

                end

                # storing the optimal value of the obejective and variables of the p-LR with fixed multipliers
                status = optimize!(bnb_node.dual_subproblems[s])

                # check whether one of the subproblems was INFEASIBLE and break the function if it was

                if termination_status_error_check(Int(termination_status(bnb_node.dual_subproblems[s])))
                    return false
                end

                dual_subrpoblems_objective_values_vector[s] = objective_value(bnb_node.dual_subproblems[s])
                fs_dv_values[ :, s ] = value.(bnb_node.dual_subproblems[s][:x])
                ss_dv_values[ :, s ] = value.(bnb_node.dual_subproblems[s][:y])
                z_values[ :, s] = value.(bnb_node.dual_subproblems[s][:z])

                if bnb_node.initial_parameters.RNMDT_is_used # if we apply LR to the RNMDT-based relaxation we store auxiliary variable w
                    RNMDT_w_values[ :, :, s] = value.(bnb_node.dual_subproblems[s][:w_RNMDT])
                end
        end # for loop

    else
        for s in 1 : bnb_node.initial_parameters.num_scen # if we don't use scenario-based parallelisation

            # objective_update

            if bnb_node.initial_parameters.RNMDT_is_used  # if we apply LR to the RNMDT-based relaxation

                if bnb_node.initial_parameters.al_is_used
                    @objective( bnb_node.dual_subproblems[s], Min,
                    - bnb_node.initial_parameters.scen_prob[s] *
                        (
                        sum( bnb_node.generated_parameters.objective_Qs[s][i, j] * bnb_node.dual_subproblems[s][:w_RNMDT][i,j]
                            for i = 1 : bnb_node.initial_parameters.num_second_stage_var,
                                j = 1 : bnb_node.initial_parameters.num_second_stage_var)
                        + sum( bnb_node.dual_subproblems[s][:x][i] * bnb_node.generated_parameters.objective_c[i]  for i = 1:bnb_node.initial_parameters.num_first_stage_var)
                        + sum( bnb_node.dual_subproblems[s][:y][j] * bnb_node.generated_parameters.objective_fs[s][j]  for j = 1:bnb_node.initial_parameters.num_second_stage_var)
                        - bnb_node.initial_parameters.μ * sum(bnb_node.dual_subproblems[s][:z][r] for r  = 1:bnb_node.initial_parameters.num_const )
                        )
                        +  sum( lm_values[s] .* bnb_node.dual_subproblems[s][:x] )
                    )


                else
                    @objective( bnb_node.dual_subproblems[s], Min,
                    - bnb_node.initial_parameters.scen_prob[s] *
                        (
                        sum( bnb_node.generated_parameters.objective_Qs[s][i, j] * bnb_node.dual_subproblems[s][:w_RNMDT][i,j]
                            for i = 1 : bnb_node.initial_parameters.num_second_stage_var,
                                j = 1 : bnb_node.initial_parameters.num_second_stage_var)
                        + sum( bnb_node.dual_subproblems[s][:x][i] * bnb_node.generated_parameters.objective_c[i]  for i = 1:bnb_node.initial_parameters.num_first_stage_var)
                        + sum( bnb_node.dual_subproblems[s][:y][j] * bnb_node.generated_parameters.objective_fs[s][j]  for j = 1:bnb_node.initial_parameters.num_second_stage_var)
                        - bnb_node.initial_parameters.μ * sum(bnb_node.dual_subproblems[s][:z][r] for r  = 1:bnb_node.initial_parameters.num_const )
                        )
                        +  sum( f_lambda_lagrangian( lm_values[:], s ) .* bnb_node.dual_subproblems[s][:x] )
                    )
                end

            else # if we apply LR straight to the primal problem

                @objective( bnb_node.dual_subproblems[s], Min,
                  - bnb_node.initial_parameters.scen_prob[s] *
                    (
                    sum( bnb_node.dual_subproblems[s][:y][i] * bnb_node.generated_parameters.objective_Qs[s][i, j] * bnb_node.dual_subproblems[s][:y][j] for i = 1 : bnb_node.initial_parameters.num_second_stage_var, j = 1 : bnb_node.initial_parameters.num_second_stage_var)
                    + sum( bnb_node.dual_subproblems[s][:x][i] * bnb_node.generated_parameters.objective_c[i]  for i = 1:bnb_node.initial_parameters.num_first_stage_var)
                    + sum( bnb_node.dual_subproblems[s][:y][j] * bnb_node.generated_parameters.objective_fs[s][j] for j = 1:bnb_node.initial_parameters.num_second_stage_var)
                    - bnb_node.initial_parameters.μ * sum(bnb_node.dual_subproblems[s][:z][r] for r  = 1:bnb_node.initial_parameters.num_const )
                    )
                    +  sum( f_lambda_lagrangian( lm_values[:], s ) .* bnb_node.dual_subproblems[s][:x] )
                )

            end

            # storing the optimal value of the obejective and variables of the p-LR with fixed multipliers
            status = optimize!(bnb_node.dual_subproblems[s])

            @show bnb_node.initial_parameters.μ * sum(value.(bnb_node.dual_subproblems[s][:z])[r] for r  = 1:bnb_node.initial_parameters.num_const )

            # check whether one of the subproblems was INFEASIBLE and break the function if it was

            if termination_status_error_check(Int(termination_status(bnb_node.dual_subproblems[s])))
                return false
            end

            dual_subrpoblems_objective_values_vector[s] = objective_value(bnb_node.dual_subproblems[s])
            fs_dv_values[ :, s ] = value.(bnb_node.dual_subproblems[s][:x])
            ss_dv_values[ :, s ] = value.(bnb_node.dual_subproblems[s][:y])
            z_values[ :, s] = value.(bnb_node.dual_subproblems[s][:z])

            if bnb_node.initial_parameters.RNMDT_is_used # if we apply LR to the RNMDT-based relaxation we store auxiliary variable w
                RNMDT_w_values[ :, :, s] = value.(bnb_node.dual_subproblems[s][:w_RNMDT])
            end

        end # for loop

    end # if (parallelisation based)
    if bnb_node.initial_parameters.al_is_used
        return [dual_subrpoblems_objective_values_vector, fs_dv_values, ss_dv_values, RNMDT_w_values, z_values]
    elseif bnb_node.initial_parameters.RNMDT_is_used # if we apply LR to the RNMDT-based relaxation we store auxiliary variable w
        return [dual_subrpoblems_objective_values_vector, fs_dv_values, ss_dv_values, RNMDT_w_values]
    else
        return [dual_subrpoblems_objective_values_vector, fs_dv_values, ss_dv_values]
    end
end

# the function checks the wether the termination status
# had any errors and returns true if it has and fals otherwise
function termination_status_error_check(status_code::Int64)
    # if INFEASIBLE or INFEASIBLE_OR_UNBOUNDED
    return  ((status_code == 2) || (status_code == 6))
end

# the function contains the implementation of the proximal bundle method
function proximal_bundle_method(bnb_node::node, zero_iteration_cg, tolerance)

## simplifying the references
    input_parameters = bnb_node.initial_parameters.bm_parameters

## initialising variables for the main loop of the proximal bundle method

    # lagrangian relaxation variables for the x and y non anticipativity conditions written in the column, for each iteration
    lm_values = Array{Array{Float64}}(undef, input_parameters.max_number_of_iterations, bnb_node.initial_parameters.num_scen - 1)
    [ lm_values[1, i] = zero_iteration_cg[i] for i = 1 : bnb_node.initial_parameters.num_scen - 1 ]

    # defining the array for storing the optimal objective function values
    # for the piecewise linear apprixmation model at each iteration
    pla_objective = Array{Float64}(undef, 1, input_parameters.max_number_of_iterations)

    # vector for storing dual objective function values calculated
    # using centre of gravity at each iteration
    do_value_cg = Array{Float64}(undef, 1, input_parameters.max_number_of_iterations)

    # subgradient vector for dual problem at each iteration
    sv_dp  = Array{Array{Float64}}(undef, input_parameters.max_number_of_iterations, bnb_node.initial_parameters.num_scen - 1)

    # the centre of gravity at each iteration
    cg = Array{Array{Float64}}(undef, input_parameters.max_number_of_iterations, bnb_node.initial_parameters.num_scen - 1)

    # vector for storing dual objective function values calculated
    # using fixed lagragian multiplers at each iteration
    do_value_lm = Array{Float64}(undef, input_parameters.max_number_of_iterations)

    # vector that contains the values of the first stage decision variables (x  in this case)
    # calucalted at each iteration
    fs_dv_values = Array{Float64}(undef, input_parameters.max_number_of_iterations, bnb_node.initial_parameters.num_first_stage_var, bnb_node.initial_parameters.num_scen)

    # vector that contains the values of the second stage decision variables (y in this case)
    # calucalted at each iteration
    ss_dv_values = Array{Float64}(undef, input_parameters.max_number_of_iterations, bnb_node.initial_parameters.num_second_stage_var, bnb_node.initial_parameters.num_scen)

    # if we cobined LR and RNMDT we create the vector that contains
    # the values of the axuliary variable w calucalted at each iteration
    if  bnb_node.initial_parameters.RNMDT_is_used
        RNMDT_w_values =  Array{Float64}(undef, input_parameters.max_number_of_iterations, bnb_node.initial_parameters.num_second_stage_var, bnb_node.initial_parameters.num_second_stage_var, bnb_node.initial_parameters.num_scen)
    end

    # subgradient vector for dual problem at each iteration
    sv_dp  = Array{Array{Float64}}(undef, input_parameters.max_number_of_iterations, bnb_node.initial_parameters.num_scen - 1)

    # initialising the master problem
    master_problem = Model(optimizer_with_attributes(Gurobi.Optimizer,  "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol, "Method" => bnb_node.initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag,  "Threads" => bnb_node.initial_parameters.gurobi_parameters.Threads, "NumericFocus" => initial_parameters.gurobi_parameters.NumericFocus, "Presolve" => 0))
    @variables master_problem begin
        θ
        lagrangian_multipliers_representing_variable[ 1 : bnb_node.initial_parameters.num_first_stage_var,
            1 : bnb_node.initial_parameters.num_scen - 1]
    end

## initialising auxiliary variables for step-size modification

    v_hat = Array{Float64}(undef, input_parameters.max_number_of_iterations)
    v = Array{Float64}(undef, input_parameters.max_number_of_iterations)
    δ_hat = Array{Float64}(undef, input_parameters.max_number_of_iterations)
    δ = Array{Float64}(undef, input_parameters.max_number_of_iterations)
    h = Array{Float64}(undef, input_parameters.max_number_of_iterations)

    # subgradient vector for master problem at each iteration
    sv_mp = Array{Array{Float64}}(undef, input_parameters.max_number_of_iterations, bnb_node.initial_parameters.num_scen - 1)

## 0 ITERATION

    # initializing the values of the lagrangian multipliers and centre of gravity for 0 iteration
    zero_iteration_lm_values = zero_iteration_cg

    # initialising the array that is going to contain the value of the stepsize
    u = Array{Float64}(undef, input_parameters.max_number_of_iterations)
    # initialising the value for the step size at the 1 iteration
    Random.seed!(bnb_node.initial_parameters.random_seed)
    u[1] =  rand(input_parameters.min_ssv : input_parameters.max_ssv)
    #@show u[1]
    # initialising the array that is going to conytain the value of the stepsize
    i = Array{Float64}(undef,  input_parameters.max_number_of_iterations)
    # initialising the value for i at the 1 iteration
    i[1] = 0

    # setting the iteration counter
    k = 0

    # solving dual subproblems using the values of the lagrangian multipliers at 0 iteration
    zero_iteration_dual_subproblems_output = dual_subproblems_values(bnb_node, zero_iteration_lm_values)

    #print(bnb_node.dual_subproblems[1])
    # check whether one of the subproblems was INFEASIBLE and break the function if it was
    #@show zero_iteration_dual_subproblems_output
    if zero_iteration_dual_subproblems_output == false
        return false
    end

    zero_iteration_do_value_cg = sum(zero_iteration_dual_subproblems_output[1])

    # calcualting the subgradient vector at 0 iteration
    zero_iteration_sv_dp = Array{Array{Float64}}(undef, bnb_node.initial_parameters.num_scen - 1 )
    [zero_iteration_sv_dp[s - 1]  =  + (zero_iteration_dual_subproblems_output[2][ :, 1] - zero_iteration_dual_subproblems_output[2][ :, s]) for  s in 2 : bnb_node.initial_parameters.num_scen ]

    # adding cuts for the master subprblem
    @constraint(master_problem, master_problem[:θ] <= sum(zero_iteration_dual_subproblems_output[1]) + sum( sum( zero_iteration_sv_dp[s] .* ( master_problem[:lagrangian_multipliers_representing_variable][:, s] .- zero_iteration_lm_values[s] ) ) for s = 1 : bnb_node.initial_parameters.num_scen - 1) )
    @objective(master_problem, Max, master_problem[:θ] - u[k+1]/2 * sum( sum( (master_problem[:lagrangian_multipliers_representing_variable][:, s] .- zero_iteration_cg[s] ).^2 ) for  s in 1 : bnb_node.initial_parameters.num_scen - 1 ) )


## MAIN LOOP
    while (k < input_parameters.max_number_of_iterations-1) && ((k == 0) ? true : v[k] > tolerance)

        # update iteration counter
        k += 1
        #@show k

        # solve the master problem and
        # gather the new values for lagrangian multipliers
        # and piecewise apprixamtion model
        optimize!(master_problem)
        [ lm_values[k, s] = value.(master_problem[:lagrangian_multipliers_representing_variable][:, s]) for s = 1 : bnb_node.initial_parameters.num_scen - 1 ]
        pla_objective[k] =  value.(master_problem[:θ])

        # calculate the subgradeint for master problem
        if k == 1
            [sv_mp[k, s - 1]  =  u[k]* (lm_values[k, s - 1] .- zero_iteration_lm_values[s - 1]) for  s in 2 : bnb_node.initial_parameters.num_scen ]
        else
            [sv_mp[k, s - 1]  =  u[k]* (lm_values[k, s - 1] .- cg[k - 1, s - 1]) for  s in 2 : bnb_node.initial_parameters.num_scen ]
        end

        # solving dual subproblems using the values of the lagrangian multipliers
        # at current iteration
        dual_subproblems_output = dual_subproblems_values(bnb_node, lm_values[k,:])
        #print(bnb_node.dual_subproblems[1])

        # check whether one of the subproblems was INFEASIBLE and break the function if it was
        if dual_subproblems_output == false
            return false
        end

        # storing dual objective function value
        # using fixed lagragian multiplers at the current iteration
        do_value_lm[k] = sum(dual_subproblems_output[1])

        print("iteraton = $k, dual value = $(do_value_lm[k])\n")
        # storing first stage decsion variables values at the current iteration
        fs_dv_values[k, :, :] = dual_subproblems_output[2]

        # storing second stage decsion variables values at the current iteration
        ss_dv_values[k, :, :] = dual_subproblems_output[3]

        # if we combined LR and RNMDT
        # storing w variables values at the current iteration
        if  bnb_node.initial_parameters.RNMDT_is_used
            RNMDT_w_values[k,:,:,:] = dual_subproblems_output[4]
        end

        # calcualting the subgradient vector at current iteration
        [sv_dp[k, s - 1]  =  + (dual_subproblems_output[2][ :, 1] - dual_subproblems_output[2][ :, s]) for  s in 2 : bnb_node.initial_parameters.num_scen]

        if k == 1
            v_hat[k] =  do_value_lm[k] - zero_iteration_do_value_cg
            v[k] =  pla_objective[k] - zero_iteration_do_value_cg
            δ_hat[k] = v_hat[k] + sum( [ dot(sv_dp[k, s], lm_values[k, s] .- zero_iteration_cg[s] ) for s = 1: bnb_node.initial_parameters.num_scen - 1 ])
            δ[k] = v[k] + sum( [ dot(sv_mp[k, s], (lm_values[k, s] .- zero_iteration_cg[s] ) ) for s = 1: bnb_node.initial_parameters.num_scen - 1 ])
            h[k] =  u[k]*(1 - ( (v[k]!=0) ? (v_hat[k]/v[k]) : 0))

        else
            v_hat[k] =  do_value_lm[k] - do_value_cg[k-1]
            v[k] =  pla_objective[k] - do_value_cg[k-1]
            δ_hat[k] = v_hat[k] + sum( [ dot(sv_dp[k, s], (lm_values[k, s] .- cg[k-1, s] ) ) for s = 1: bnb_node.initial_parameters.num_scen - 1 ])
            δ[k] = v[k] + sum( [ dot(sv_mp[k, s], (lm_values[k, s] .- cg[k-1, s] ) ) for s = 1: bnb_node.initial_parameters.num_scen - 1 ])
            h[k] =  u[k]*(1 - ( (v[k]!=0) ? (v_hat[k]/v[k]) : 0) )
        end
        #@show v_hat[k]
        #@show v[k]
        if v_hat[k] >= input_parameters.mL*v[k]
            cg[k, :] = lm_values[k, :]
            do_value_cg[k] = do_value_lm[k]
            #print("SERIOUS STEP\n")

            # STEPSIZE UPDATE
            u[k+1] = u[k]
            #@show h[k]
            if (v_hat[k] >= input_parameters.mR * v[k]) && (i[k]>0)
                u[k+1] = max(h[k], 0.1 * u[k], input_parameters.min_ssv)
            elseif (i[k]>3)
                u[k+1] = max( 0.5 * u[k], input_parameters.min_ssv)
            end

            if (u[k+1] == u[k])
                i[k+1] = max(i[k]+1, 1)
            else
                i[k+1] = 1
            end

        else
            if k == 1
                cg[k, :] = zero_iteration_cg
                do_value_cg[k] = zero_iteration_do_value_cg
            else
                cg[k, :] = cg[k-1, :]
                do_value_cg[k] = do_value_cg[k-1]
            end

            # STEPSIZE UPDATE
            u[k+1] = u[k]
            if (δ_hat[k] > max(δ[k] + sqrt( sum( dot(sv_mp[k,s],sv_mp[k,s]) for s = 1 : bnb_node.initial_parameters.num_scen - 1 ) ), 10 * v[k]) ) && (i[k] < -3)
                u[k+1] = min(h[k], 10 * u[k])
            else
                u[k+1] = min(10 * u[k], input_parameters.max_ssv)
            end

            if (u[k+1] == u[k])
                i[k+1] = min(i[k] - 1, -1)
            else
                i[k+1] = -1
            end

        end
    @show u[k+1]
    #@show i[k]

    @constraint(master_problem, master_problem[:θ] <= sum(dual_subproblems_output[1]) + sum( sum( sv_dp[k, s] .* ( master_problem[:lagrangian_multipliers_representing_variable][:, s] .- lm_values[k, s] ) ) for s = 1 : bnb_node.initial_parameters.num_scen - 1) )
    @objective(master_problem, Max, master_problem[:θ] - u[k+1]/2 * sum( sum( (master_problem[:lagrangian_multipliers_representing_variable][:, s] .- cg[k, s] ).^2 ) for  s in 1 : bnb_node.initial_parameters.num_scen - 1 ) )
    #print(master_problem)
    end

    if bnb_node.initial_parameters.RNMDT_is_used # if we apply LR to the RNMDT-based relaxation we store auxiliary variable w
        return dm_output(do_value_cg[1:k], [fs_dv_values[k, :, :], ss_dv_values[k, :, :], RNMDT_w_values[k, :, :, :]])
    else
        return dm_output(do_value_cg[1:k], [fs_dv_values[k, :, :], ss_dv_values[k, :, :]])
    end

end



##
