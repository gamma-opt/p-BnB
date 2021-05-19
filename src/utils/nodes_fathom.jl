"""
    nodes_fathom(X::bnb_model)
If the first candidate for the optimal solutin is found
fathom the nodes whose parent dual bounds are higher than current optimal solution

"""
function nodes_fathom(X::bnb_model)

    # find the positions of the nodes in the stack
    # for which the pn_db > UBDg
    fathom_indices = findall(x -> x > X.UBDg, X.pn_db)

    #print("node_fathomed\n")
    #@show X.UBDg
    #@show X.pn_db[fathom_indices]

    # fathoming the nodes at correspondent position
    splice!(X.nodes, fathom_indices)

    # fathoming the correspondent parent dual values
    splice!(X.pn_db, fathom_indices)


end
