"""
    a node_copy(parent_node::node)
is an auxiliary function to copy all the fileds of parent node into the child node

"""

function node_copy(parent_node::node)

    #simplifying the notations
    initial_parameters = parent_node.initial_parameters

    child_node = node(copy(parent_node.primal_problem), [copy(parent_node.dual_subproblems[j]) for j = 1:length(parent_node.dual_subproblems)], parent_node.initial_parameters, parent_node.generated_parameters)

    @suppress set_optimizer(child_node.primal_problem, optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol,
    "Method" => initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag, "Threads" => initial_parameters.gurobi_parameters.Threads))

    @suppress [set_optimizer(child_node.dual_subproblems[j], optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol,
    "Method" => initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag, "Threads" => initial_parameters.gurobi_parameters.Threads))
             for j = 1:length(child_node.dual_subproblems)]

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
    
    if initial_parameters.pool_problem_is_used
        
        model, bub = pooling_problem_one_layer_pools(initial_parameters.pool_prob_par.primal_pool_problem_link, initial_parameters.pool_prob_par.start_value_pooling_nodes)
        
        initial_parameters, generated_parameters, var, bub = polling_problem_parameters_generation(model, initial_parameters.num_scen, initial_parameters.pool_prob_par.flow_cost, initial_parameters.pool_prob_par.pool_cost, initial_parameters.RNMDT_precision_factor[1], bub)
        primal_problem = pooling_MIP_generation(initial_parameters, generated_parameters, var, bub, initial_parameters.pool_prob_par.stoc_pp)

        dual_subproblems = RNMDT_based_augmented_lagrangian_relaxation_pooling_problem_generation(initial_parameters, generated_parameters, var, bub, initial_parameters.pool_prob_par.stoc_pp)

    else 
    
        generated_parameters = parameters_generation(initial_parameters)
        primal_problem = MIP_generation(initial_parameters, generated_parameters)

        if initial_parameters.RNMDT_is_used
            if initial_parameters.al_is_used
                dual_subproblems = RNMDT_based_augmented_lagrangian_relaxation_problem_generation(initial_parameters, generated_parameters, initial_parameters.RNMDT_precision_factor)
            else
                dual_subproblems = RNMDT_based_lagrangian_relaxation_problem_generation(initial_parameters, generated_parameters, initial_parameters.RNMDT_precision_factor)
            end
        else
            dual_subproblems = primal_problem_based_lagrangian_relaxation_generation(initial_parameters, generated_parameters)
        end
    end
    #@show initial_parameters.num_first_stage_var
    zero_node = node(primal_problem, dual_subproblems, initial_parameters, generated_parameters)
        
    return zero_node, initial_parameters
end
