"""
    node
Stores structers used in stack nodes. Has the following fields:
* `primal_problem::Model`:                      primal problem associated with the node
* `dual_subproblems::Vector{Model}`:            vector of dual subpeoblems associated with the node
"""
mutable struct node
    primal_problem::Model
    dual_subproblems::Vector{Model}
end
