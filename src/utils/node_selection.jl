"""
    node_selection(X::bnb_model)
Selects node with the biggest upper bound. Returns `(Node,UBD,id)`
where `Node` is the correspodent node, `UBD` is the
upper bound of the node, and `id` is the id number of the node.

"""
function node_selection(X::bnb_model)
    return splice!(X.Nodes,length(X.Nodes)), splice!(X.id,length(X.id))
end
