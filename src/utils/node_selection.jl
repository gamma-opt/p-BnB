"""
    node_selection(X::bnb_model)
Selects node with the biggest upper bound. Returns `(Node,UBD,id)`
where `Node` is the correspodent node, `UBD` is the
upper bound of the node, and `id` is the id number of the node.

"""
function node_selection(X::bnb_model)

    # first, find the minimum among all the parend dual bond values associated with each node.
    min_pn_db = minimum(X.pn_db)

    # find the frist position from the top of the stack of the node the pn_db  =  min_pn_db
    # (it is done so that if we happen to have a few nodes with pn_db  =  min_pn_db than we will go with
    # the most right-hand sided -> the one we added to the stack the lattest since the
    # right-hand side branching problem is always added to the stack after left-hand side branching)

    selcted_id = findlast( x -> x == min_pn_db, X.pn_db)

    print("node selected with pn_db $(X.pn_db[selcted_id])")
    @show X.pn_db

    return splice!(X.nodes, selcted_id), splice!(X.pn_db, selcted_id)
    #return splice!(X.nodes,length(X.nodes)), splice!(X.id,length(X.id))
end
