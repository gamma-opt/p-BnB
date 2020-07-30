"""
    bnb_solve(intial_parameters::MIP_initial_parameters, epsilon::Float64, epsilon_bb::Float64)
returns the resulting bnb_struct of caroe and shultz bnb method with
the parameters epsilon and epsilon_bb. The method is applied to the optimisation
problem formulated using intial_parameters.
"""

function bnb_solve(intial_parameters::MIP_initial_parameters, epsilon::Float64, epsilon_bb::Float64)

    # to use Ipopt for the primal problems
    # when calculating LB
    intial_parameters.is_int_fixed = false

    # define the starting centre of gravity for the bundle method
    initial_center_of_gravity = Array{Any}(undef, initial_parameters.num_scen - 1)
    [ initial_center_of_gravity[i] = zeros(1,initial_parameters.num_int_var)
        for i = 1 : initial_parameters.num_scen - 1 ]

    init_time = time()

    # create first bnb_structure
    bnb_struct = bnb_model(initial_parameters)

    while length(bnb_struct.nodes)>0

        #2 extract the node from the stack and try to optimize it
        current_node, current_node_id = node_selection(bnb_struct)

        # using the budnle method to calculate the dual value of the chosen node
        bm_output = bundle_method(current_node, 100, initial_center_of_gravity, "parallelised" )

        # compute aaveraged values of theinteger variables
        current_avg_int_var = sum(current_node.initial_parameters.scen_prob[j] .* bm_output.int_var[:,j] for j = 1:current_node.initial_parameters.num_scen)

        # calcuate the delta
        current_delta = delta_computation(bm_output.int_var)

        # if the values of int var coincide for all the scenrios
        if maximum(current_delta) == 0

            # save new dual value if it is lower than exisiting
            if bm_output.dual_objective_value[end]<bnb_struct.UBDg

                bnb_struct.UBDg = bm_output.dual_objective_value[end]
                push!(bnb_struct.UBDg_hist, bnb_struct.UBDg)
                push!(bnb_struct.UBDg_time_hist, time()-init_time)

                bnb_struct.max_id = current_node_id

                bnb_struct.first_fnd = true
                bnb_struct.soln_val = current_avg_int_var
                push!(bnb_struct.soln_hist, current_avg_int_var)

            end

        elseif maximum(current_delta) > epsilon && bm_output.dual_objective_value[end]<bnb_struct.UBDg

            # if some of the integer variables are still continuous
            if sum(isinteger.(current_avg_int_var))/length(current_avg_int_var) == 1

                # find the first index that breaks integrality conditions
                variable_index = findfirst(isequal(0), isinteger.(current_avg_int_var))

                # do branching on that index to ensure integrality conditions
                l_node, r_node = branching(current_node, current_avg_int_var, variable_index, "int")

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
                l_node, r_node = branching(current_node, current_avg_int_var, variable_index, "cont", epsilon_bb)

                push!(bnb_struct.nodes, l_node)
                push!(bnb_struct.id, length(bnb_struct.nodes)+1)

                push!(bnb_struct.nodes, r_node)
                push!(bnb_struct.id, length(bnb_struct.nodes)+1)

            end #if

        end #if

        # update LB

        new_LB = update_LB(current_node, current_avg_int_var)
        if new_LB > bnb_struct.LBDg
            bnb_struct.LBDg = new_LB
            push!(bnb_struct.LBDg_hist, new_LB)
            push!(bnb_struct.LBDg_time_hist, time() - init_time)
        end

    end #while

    return bnb_struct

end #function
