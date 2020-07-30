function update_LB(current_node::node, int_var_avg::Vector{Float64})

 for i = 1:current_node.initial_parameters.num_scen
     fix.(current_node.primal_problem[:x][:,i], int_var_avg, force = true)
 end

 optimize!(current_node.primal_problem)

 return(objective_value(current_node.primal_problem))

end
