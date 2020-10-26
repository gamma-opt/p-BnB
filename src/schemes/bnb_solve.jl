"""
    bnb_solve(initial_parameters::MIP_initial_parameters, non_ant_tol::Float64, tol_bb::Float64, integrality_tolerance::Float64)
returns the resulting bnb_struct of caroe and shultz bnb method with
the parameters non_ant_tol and tol_bb. The method is applied to the optimisation
problem formulated using initial_parameters.
"""

function bnb_solve(initial_parameters::MIP_initial_parameters, non_ant_tol::Float64, tol_bb::Float64, integrality_tolerance::Float64)

    # define the starting centre of gravity for the bundle method
    initial_centre_of_gravity = intial_centre_of_gravity_generation(initial_parameters.bm_parameters.initial_centre_of_gravity_min, initial_parameters.bm_parameters.initial_centre_of_gravity_max, initial_parameters.random_seed, initial_parameters.num_scen, initial_parameters.num_first_stage_var)

    init_time = time()

    # create first bnb_structure
    bnb_struct = bnb_model(initial_parameters)

    while length(bnb_struct.nodes)>0

        print("length before = $(length(bnb_struct.nodes))\n")

        #2 extract the node from the stack and try to optimize it
        current_node, current_node_id = node_selection(bnb_struct)

        print("length after = $(length(bnb_struct.nodes))\n")

        # using the budnle method to calculate the dual value of the chosen node
        time_init = time()
        # if initial_parameters.RNMDT_is_used
        #    current_node_ouput = dynamic_precision_RNMDT_algorithm(current_node)
        # else
        current_node_ouput = bundle_method(current_node, initial_centre_of_gravity)
        # end
        print("node id: $current_node_id, Dual bound: $(current_node_ouput.dual_objective_value[end])\n")
        print("node id: $current_node_id, first stage variables: $(current_node_ouput.variables_values[1])\n")
        # compute averaged values of the first stage variables
        current_avg_first_stage_var = round.(sum(current_node.initial_parameters.scen_prob[j] .* current_node_ouput.variables_values[1][:,j] for j = 1:current_node.initial_parameters.num_scen), digits = 5)

        # calcuate the delta
        current_delta = delta_computation(current_node_ouput.variables_values[1])

        # if the values of integer variables coincide for all the scenrios (implying the solution is feasible for the primal problem)
        if maximum(current_delta) == 0

            # save new dual value as the global primal (lower) bound if it is higher than exisiting
            if current_node_ouput.dual_objective_value[end]>bnb_struct.LBDg

                bnb_struct.LBDg = current_node_ouput.dual_objective_value[end]
                push!(bnb_struct.LBDg_hist, bnb_struct.LBDg)
                push!(bnb_struct.LBDg_time_hist, time()-init_time)

                bnb_struct.max_id = current_node_id

                bnb_struct.first_fnd = true
                bnb_struct.soln_val = current_avg_first_stage_var
                push!(bnb_struct.soln_hist, current_avg_first_stage_var)

                # update global UB if the node was not fathomed (if current_node_ouput.dual_objective_value[end]>bnb_struct.LBDg )
                # and the dual bound at the current iteration is smaller than exsiting global upper bound
                if  current_node_ouput.dual_objective_value[end] < bnb_struct.UBDg
                    bnb_struct.UBDg = current_node_ouput.dual_objective_value[end]
                    push!(bnb_struct.UBDg_hist, bnb_struct.UBDg)
                end

            end

        elseif maximum(current_delta) > non_ant_tol && current_node_ouput.dual_objective_value[end]>bnb_struct.LBDg

            # if some of the first stage variables with integer idexes are still continuous
            if sum(isinteger.(current_avg_first_stage_var[current_node.generated_parameters.x_int_indexes]))/length(current_avg_first_stage_var[current_node.generated_parameters.x_int_indexes]) != 1

                # find the first index that breaks integrality conditions
                variable_index = current_node.generated_parameters.x_int_indexes[findfirst(isequal(0), isinteger.(current_avg_first_stage_var[current_node.generated_parameters.x_int_indexes]))]

                print("current $(current_avg_first_stage_var)\n")
                print("variable index is  $variable_index, value is $(current_avg_first_stage_var[variable_index])\n")

                # do branching on that index to ensure integrality conditions
                l_node, r_node = branching(current_node, current_avg_first_stage_var, variable_index, "int")

                push!(bnb_struct.nodes, l_node)
                push!(bnb_struct.id, length(bnb_struct.nodes)+1)

                push!(bnb_struct.nodes, r_node)
                push!(bnb_struct.id, length(bnb_struct.nodes)+1)

            else

                # otherwise, find the index that breaks non-anticipativy
                # conditions the most
                variable_index = argmax(current_delta)

                # do branching on that index to ensure
                # non-anticipativity conditions
                l_node, r_node = branching(current_node, current_avg_first_stage_var, variable_index, "cont", tol_bb)

                push!(bnb_struct.nodes, l_node)
                push!(bnb_struct.id, length(bnb_struct.nodes)+1)

                push!(bnb_struct.nodes, r_node)
                push!(bnb_struct.id, length(bnb_struct.nodes)+1)

            end #if

            # update global UB if the node was not fathomed (if current_node_ouput.dual_objective_value[end]>bnb_struct.LBDg )
            # and the dual bound at the current iteration is smaller than exsiting global upper bound
            if current_node_ouput.dual_objective_value[end] < bnb_struct.UBDg
                bnb_struct.UBDg = current_node_ouput.dual_objective_value[end]
                push!(bnb_struct.UBDg_hist, bnb_struct.UBDg)
            end


        end #elseif


    end # while

    return bnb_struct

end #function
