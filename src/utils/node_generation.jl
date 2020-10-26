"""
    a node_copy(parent_node::node)
is an auxiliary function to copy all the fileds of parent node into the child node

"""

function node_copy(parent_node::node)
    child_node = node(copy(parent_node.primal_problem), [copy(parent_node.dual_subproblems[j]) for j = 1:length(parent_node.dual_subproblems)], parent_node.initial_parameters, parent_node.generated_parameters)
    @suppress set_optimizer(child_node.primal_problem, optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2))
    @suppress [set_optimizer(child_node.dual_subproblems[j], optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2)) for j = 1:length(child_node.dual_subproblems)]

    return child_node
end

"""
    a child_node_generation(parent_node::node, var_index::Int, inequality::String, value::Int)
returns a child node generated from a parent_node by copying it
and addinbg auxiliary constraint x[var_index] "<=" or ">=" value
to primal problem and dual subproblems.

"""
function child_node_generation(parent_node::node, var_index::Int, inequality::String, value::Float64)
    child_node = node_copy(parent_node)
    #child_node = node(copy(parent_node.primal_problem), [copy(parent_node.dual_subproblems[j]) for j = 1:length(parent_node.dual_subproblems)], parent_node.initial_parameters, parent_node.generated_parameters)
    #@suppress set_optimizer(child_node.primal_problem, optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2))
    inequality == "<=" ? @constraint(child_node.primal_problem, child_node.primal_problem[:x][var_index, :] .<= value) : @constraint(child_node.primal_problem, child_node.primal_problem[:x][var_index, :] .>= value)
    #@suppress [set_optimizer(child_node.dual_subproblems[j], optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2)) for j = 1:length(child_node.dual_subproblems)]
    [ inequality == "<=" ? @constraint(child_node.dual_subproblems[j], child_node.dual_subproblems[j][:x][var_index] <= value) : @constraint(child_node.dual_subproblems[j], child_node.dual_subproblems[j][:x][var_index] >= value) for j = 1:length(child_node.dual_subproblems)]
    return child_node
end

"""
    zero_node_generation(initial_parameters::MIP_initial_parameters)
returns zero node by relaxing integrality conditions
in the original problem

"""
function zero_node_generation(initial_parameters::MIP_initial_parameters)
    generated_parameters = parameters_generation(initial_parameters)
    primal_problem = MIP_generation(initial_parameters, generated_parameters)
    if initial_parameters.RNMDT_is_used
        dual_subproblems = RNMDT_based_lagrangian_relaxation_problem_generation(initial_parameters, generated_parameters, initial_parameters.RNMDT_precision_factor)
    else
        dual_subproblems = primal_problem_based_lagrangian_relaxation_generation(initial_parameters, generated_parameters)
    end
    zero_node = node(primal_problem, dual_subproblems, initial_parameters, generated_parameters)
    return zero_node
end
