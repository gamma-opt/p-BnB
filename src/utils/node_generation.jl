"""
    a child_node_generation(parent_node::node, var_index::Int, inequality::String, value::Int)
returns a child node generated from a parent_node by copying it
and addinbg auxiliary constraint x[var_index] "<=" or ">=" value
to primal problem and dual subproblems.

"""
function child_node_generation(parent_node::node, var_index::Int, inequality::String, value::Float64)
    child_node = node(copy(parent_node.primal_problem), [copy(parent_node.dual_subproblems[j]) for j = 1:length(parent_node.dual_subproblems)], parent_node.initial_parameters, parent_node.generated_parameters)
    @suppress set_optimizer(child_node.primal_problem, optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2))
    inequality == "<=" ? @constraint(child_node.primal_problem, child_node.primal_problem[:x][var_index, :] .<= value) : @constraint(child_node.primal_problem, child_node.primal_problem[:x][var_index, :] .>= value)
    @suppress [set_optimizer(child_node.dual_subproblems[j], optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2)) for j = 1:length(child_node.dual_subproblems)]
    [ inequality == "<=" ? @constraint(child_node.dual_subproblems[j], child_node.dual_subproblems[j][:x][var_index] <= value) : @constraint(child_node.dual_subproblems[j], child_node.dual_subproblems[j][:x][var_index] >= value) for j = 1:length(child_node.dual_subproblems)]
    return child_node
end

"""
    zero_node_generation(initial_parameters::MIP_initial_parameters)
returns zero node by relaxing integrality conditions
in the original problem

"""
function zero_node_generation(initial_parameters::MIP_initial_parameters)
    primal_problem = MIP_generation(initial_parameters)
    #JuMP.unset_integer.(primal_problem[:x])
    dual_subproblems = MIP_lagrangian_relaxation_generation(initial_parameters)
    [JuMP.unset_integer.(dual_subproblems[i][:x]) for i = 1:length(dual_subproblems)]
    zero_node = node(primal_problem, dual_subproblems, initial_parameters, parameters_generation(initial_parameters))
    return zero_node
end
