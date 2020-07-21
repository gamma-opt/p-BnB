function bnb_solve(parameters::initial_parameters)

    #1   create bnb_structure
    bnb_struct = bnb_model(parameters)

    while length(bnb_struct.nodes)>0

        #2 extract the node from the stack and try to optimize it
        current_node, current_node_id = node_selection(bnb_struct)

        @suppress optimize!(current_node.)


        #3 if the solution is feasible
        if has_values(current_node) != false
            # 3.1 check the whether the node UB is not worse than the global UB
            if objective_value(current_node) >= bnb_struct.UBDg
                # 3.1.1  check the integrality conditions for each variable
                is_all_integer = true
                variable_index = 0

                while is_all_integer && variable_index < length(current_node[:x])
                    variable_index += 1
                    is_all_integer = isinteger(value(current_node[:x][variable_index]))
                end #while

                if is_all_integer
                    bnb_struct.UBDg = objective_value(current_node)
                    bnb_struct.first_fnd = true
                    bnb_struct.soln_val = value.(current_node[:x])
                else

                    l_node, r_node = branching(current_node, variable_index)

                    push!(bnb_struct.Nodes, l_node)
                    push!(bnb_struct.id, length(bnb_struct.Nodes)+1)

                    push!(bnb_struct.Nodes, r_node)
                    push!(bnb_struct.id, length(bnb_struct.Nodes)+1)
                end #if

            end #if

        end #if

    end #while

    if bnb_struct.first_fnd
        print("optmal objective value: $(bnb_struct.UBDg)\n")
        print("optimal decision variables value: $(bnb_struct.soln_val)\n")
    else print("optimal solution was not found")
    end

    return bnb_struct

end
