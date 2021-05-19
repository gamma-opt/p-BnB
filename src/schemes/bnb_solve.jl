"""
    bnb_solve(initial_parameters::MIP_initial_parameters, non_ant_tol::Float64, tol_bb::Float64, integrality_tolerance::Float64)
returns the resulting bnb_struct of caroe and shultz bnb method with
the parameters non_ant_tol and tol_bb. The method is applied to the optimisation
problem formulated using initial_parameters.
"""

function bnb_solve(initial_parameters::MIP_initial_parameters, non_ant_tol::Float64, tol_bb::Float64, integrality_tolerance::Float64)

    # define the starting centre of gravity for the bundle method
    initial_centre_of_gravity = intial_centre_of_gravity_generation(initial_parameters.bm_parameters.initial_centre_of_gravity_min, initial_parameters.bm_parameters.initial_centre_of_gravity_max, initial_parameters.random_seed, (initial_parameters.al_is_used ? initial_parameters.num_scen : initial_parameters.num_scen - 1), initial_parameters.num_first_stage_var)
    #@show initial_centre_of_gravity
    init_time = time()

    # create first bnb_structure
    bnb_struct = bnb_model(initial_parameters)
    number_of_nodes_used = 0
    #i = 1


    bnb_FW_iterations_output = []

    while length(bnb_struct.nodes) > 0
    #while i <= 1
        number_of_nodes_used += 1

        #print("\n\n number_of_nodes_used =  $number_of_nodes_used \n")
        bnb_struct.nodes_used = number_of_nodes_used
        #print("length before = $(length(bnb_struct.nodes))\n")

        #2 extract the node from the stack and try to optimize it
        current_node, current_node_id = node_selection(bnb_struct)

        #print("length after = $(length(bnb_struct.nodes))\n")

        # using the budnle method to calculate the dual value of the chosen node
        time_init = time()
        # if initial_parameters.RNMDT_is_used
        #    current_node_output = dynamic_precision_RNMDT_algorithm(current_node)
        # else


        if initial_parameters.al_is_used
            # generating intiial set V_0
            V_0, x_0 = FW_PH_V_0_initialisation(current_node, initial_centre_of_gravity)

            x_0_penalty = Array{Float64}(undef, initial_parameters.num_first_stage_var, initial_parameters.num_scen)
            [x_0_penalty[:,s] = x_0[s] for s = 1:initial_parameters.num_scen]
            current_node.initial_parameters.al_penalty_parameter = penalty_parameter_update(current_node.initial_parameters, current_node.generated_parameters, x_0_penalty)

            #@show V_0
            x_k, y_k, w_RNMDT_k, z_FR_k, z_k, w_k, ϕ_k, dual_feasibility_condition, primal_dual_residial, feasibility_set = FW_PH(current_node, V_0, x_0, initial_centre_of_gravity, number_of_nodes_used)
            #@show z_FR_k
            x_0 = Array{Float64}(undef, initial_parameters.num_first_stage_var, initial_parameters.num_scen)
            [x_0[:,s] = x_k[s] for s = 1:initial_parameters.num_scen]

            y_0 = Array{Float64}(undef, initial_parameters.num_second_stage_var, initial_parameters.num_scen)
            [y_0[:,s] = y_k[s] for s = 1:initial_parameters.num_scen]

            w_RNMDT_0 = Array{Float64}(undef, initial_parameters.num_second_stage_var, initial_parameters.num_second_stage_var, initial_parameters.num_scen)
            [w_RNMDT_0[:,:,s] = w_RNMDT_k[s] for s = 1:initial_parameters.num_scen]

            #@show x_0
            current_node_output = dm_output((ϕ_k), [x_0, y_0, w_RNMDT_0])
        else
            current_node_output = proximal_bundle_method(current_node, initial_centre_of_gravity, 1e-9)

        end

        # if we need to plot the FW iterations and we are at the root node
        if initial_parameters.al_is_used && initial_parameters.PH_SDM_parameters.PH_plot && (number_of_nodes_used==1)
            bnb_FW_iterations_output = [copy(current_node_output.dual_objective_value), dual_feasibility_condition, primal_dual_residial, feasibility_set]
        end

        #print(current_node.dual_subproblems[1])
        #print(current_node.dual_subproblems[2])
        # if one of the subproblems was infeasible, fathom the node, throw an exception and
        # move on to the next iteration immediately
        if current_node_output == false
            print("one of the subproblems was INFEASIBLE: the node is fathomed\n")
            continue
        end
        # end
        #print("node id: $current_node_id, Dual bound: $(current_node_output.dual_objective_value[end])\n")
        #print("node id: $current_node_id, first stage variables: $(current_node_output.variables_values[1])\n")
        # compute averaged values of the first stage variables
        #@show current_node.initial_parameters.num_scen
        #@show current_node.initial_parameters.scen_prob
        #@show  current_node_output.variables_values[1]

        current_avg_first_stage_var = round.(sum(current_node.initial_parameters.scen_prob[j] .* current_node_output.variables_values[1][:,j] for j = 1:current_node.initial_parameters.num_scen), digits = 8)
        #@show current_avg_first_stage_var
        #@show [current_node.initial_parameters.scen_prob[j] for j = 1:current_node.initial_parameters.num_scen]
        # calcuate the delta
        current_delta = round.(delta_computation(current_node_output.variables_values[1]), digits = 8)
        #@show current_delta
        # if the values of integer variables coincide for all the scenrios (implying the solution is feasible for the primal problem)
        #if maximum(current_delta) == 0
        if maximum(current_delta) <= non_ant_tol

            # if some of the first stage variables with integer idexes are still continuous
            if initial_parameters.al_is_used && !isempty(current_node.generated_parameters.x_int_indexes) && (sum(isinteger.(current_avg_first_stage_var[current_node.generated_parameters.x_int_indexes]))/length(current_avg_first_stage_var[current_node.generated_parameters.x_int_indexes]) != 1)

                # find the first index that breaks integrality conditions
                variable_index = current_node.generated_parameters.x_int_indexes[findfirst(isequal(0), isinteger.(current_avg_first_stage_var[current_node.generated_parameters.x_int_indexes]))]

                # do branching on that index to ensure integrality conditions
                l_node, r_node = branching(current_node, current_avg_first_stage_var, variable_index, "int")

                push!(bnb_struct.nodes, l_node)
                #push!(bnb_struct.id, length(bnb_struct.nodes)+1)
                push!(bnb_struct.pn_db, current_node_output.dual_objective_value[end])

                push!(bnb_struct.nodes, r_node)
                #push!(bnb_struct.id, length(bnb_struct.nodes)+1)
                push!(bnb_struct.pn_db, current_node_output.dual_objective_value[end])

            else # if all the first stage variables with integer idexes are integer

                # save new dual value as the global primal (upper) bound if it is lower than exisiting
                if current_node_output.dual_objective_value[end]<bnb_struct.UBDg

                    bnb_struct.UBDg = current_node_output.dual_objective_value[end]
                    push!(bnb_struct.UBDg_hist, bnb_struct.UBDg)
                    push!(bnb_struct.UBDg_time_hist, time()-init_time)

                    bnb_struct.first_fnd = true
                    bnb_struct.soln_val = current_avg_first_stage_var
                    bnb_struct.RNMDT_gap_wy = RNMDT_gap_computation(current_node_output.variables_values[2], current_node_output.variables_values[3])
                    push!(bnb_struct.soln_hist, current_avg_first_stage_var)

                    # update global LB if the node was not fathomed (if current_node_output.dual_objective_value[end]>bnb_struct.LBDg )
                    # and the dual bound at the current iteration is bigger than exsiting global lower bound
                    if  current_node_output.dual_objective_value[end] > bnb_struct.LBDg
                        bnb_struct.LBDg = current_node_output.dual_objective_value[end]
                        push!(bnb_struct.LBDg_hist, bnb_struct.LBDg)
                    end

                    # fathom the nodes who parent dual values are higher than current UBDg
                    nodes_fathom(bnb_struct)

                 end

            end

        elseif (maximum(current_delta) > non_ant_tol) && (current_node_output.dual_objective_value[end]<bnb_struct.UBDg)

            # if some of the first stage variables with integer idexes are still continuous
            if !isempty(current_node.generated_parameters.x_int_indexes) && (sum(isinteger.(current_avg_first_stage_var[current_node.generated_parameters.x_int_indexes]))/length(current_avg_first_stage_var[current_node.generated_parameters.x_int_indexes]) != 1)

                # find the first index that breaks integrality conditions
                variable_index = current_node.generated_parameters.x_int_indexes[findfirst(isequal(0), isinteger.(current_avg_first_stage_var[current_node.generated_parameters.x_int_indexes]))]

                #print("current $(current_avg_first_stage_var)\n")
                #print("variable index is  $variable_index, value is $(current_avg_first_stage_var[variable_index])\n")

                # do branching on that index to ensure integrality conditions
                l_node, r_node = branching(current_node, current_avg_first_stage_var, variable_index, "int")

                push!(bnb_struct.nodes, l_node)
                #push!(bnb_struct.id, length(bnb_struct.nodes)+1)
                push!(bnb_struct.pn_db, current_node_output.dual_objective_value[end])

                push!(bnb_struct.nodes, r_node)
                #push!(bnb_struct.id, length(bnb_struct.nodes)+1)
                push!(bnb_struct.pn_db, current_node_output.dual_objective_value[end])

            else

                # otherwise, find the index that breaks non-anticipativy
                # conditions the most
                variable_index = argmax(current_delta)

                #@show variable_index
                #@show current_avg_first_stage_var
                # do branching on that index to ensure
                # non-anticipativity conditions
                l_node, r_node = branching(current_node, current_avg_first_stage_var, variable_index, "cont", tol_bb)

                push!(bnb_struct.nodes, l_node)
                #push!(bnb_struct.id, length(bnb_struct.nodes)+1)
                push!(bnb_struct.pn_db, current_node_output.dual_objective_value[end])

                push!(bnb_struct.nodes, r_node)
                #push!(bnb_struct.id, length(bnb_struct.nodes)+1)
                push!(bnb_struct.pn_db, current_node_output.dual_objective_value[end])

            end #if

            # update global LB if the node was not fathomed (if current_node_output.dual_objective_value[end]>bnb_struct.LBDg )
            # and the dual bound at the current iteration is bigger than exsiting global lower bound
            if (current_node_output.dual_objective_value[end]<bnb_struct.UBDg) && (current_node_output.dual_objective_value[end] > bnb_struct.LBDg)
                bnb_struct.LBDg = current_node_output.dual_objective_value[end]
                push!(bnb_struct.LBDg_hist, bnb_struct.LBDg)
            end


        end #elseif


    end # while


    # if we need to plot the FW iterations and we are at the root node
    if initial_parameters.al_is_used && initial_parameters.PH_SDM_parameters.PH_plot
        return (bnb_struct, bnb_FW_iterations_output)
    else
        return bnb_struct
    end
end #function
