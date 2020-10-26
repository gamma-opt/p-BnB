"""
    node
Stores structures related to the stack nodes. Has the following fields:
* `primal_problem::Model`:                      primal problem associated with the node
* `dual_subproblems::Vector{Model}`:            vector of dual subproblems associated with the node
                                                (they can be either genereted by applying LR
                                                or RNMDT+LR to the primal problem depending
                                                on what is indicated in the structure
                                                MIP_initial_paramters filed is_RNMDT_used)

"""
mutable struct node
    primal_problem::Model
    dual_subproblems::Vector{Model}
    initial_parameters::MIP_initial_parameters
    generated_parameters::MIP_generated_parameters
end
