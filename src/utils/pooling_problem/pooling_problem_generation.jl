## This file contains all the necessary fucntions to define the pooling problem. 

function pooling_problem_one_layer_pools(src_link::String, pool_start_value::Array{Float64})
    ## Read arc data
    adata = CSV.read(src_link*"/arcs.csv", DataFrame) # Julia has base function read(); CSV.read tell we want the CSV's read(function); 
    #println(adata)

    ## Set arc array
    inode = adata[!, 1]
    jnode = adata[!, 2]
    A     = Arc.(inode, jnode)

    ## Read the node data 
    ndata = CSV.read(src_link*"/nodes.csv", DataFrame)
    #println(ndata)

    ## Number of rows in the nodes.csv file
    nnode = size(ndata, 1)   

    ## Split V into source, pool, and target nodes 
    V     = string.(ndata[!, 1])                                 # All nodes                   
    S     = [V[i] for i = 1:nnode if ndata[i,2] == "s"]  # Source nodes
    P     = [V[i] for i = 1:nnode if ndata[i,2] == "p"]  # Pool nodes
    T     = [V[i] for i = 1:nnode if ndata[i,2] == "t"]  # Target nodes

    ## Resource bounds of each node
    cost  = ndata[!, 3]
    qlb   = ndata[!, 4]    # Lower bounds of property 1
    qub   = ndata[!, 5]    # Upper bounds of property 1
    bub   = ndata[!, 6]    # Upper flow bounds

    ## Set node array and make a dictionary for convenience  
    N = Node.(cost, qlb, qub, bub)
    N = Dict(V[i] => N[i] for i = 1:size(V,1))
    #@show N

    ## Define solver
    pooling_problem = Model(() -> Gurobi.Optimizer(GRB_ENV))
    @suppress set_optimizer_attribute(pooling_problem, "NonConvex", 2)

    ## Variables
    @variable(pooling_problem, x[A] >= 0)       # Arc flows
    @variable(pooling_problem, q[V] >= 0)  # Property values

    ## Cost and revenue
    @expression(pooling_problem, cost,    sum(N[s].cost*x[a] for s in S, a in A if a.i == s))
    @expression(pooling_problem, revenue, sum(N[t].cost*x[a] for t in T, a in A if a.j == t))

    ## Objective
    @objective(pooling_problem, Max, revenue - cost)

    ## Constraints
    # Flow balance
    @constraint(pooling_problem, [p in P], sum(x[a] for a in A if a.j == p) == sum(x[a] for a in A if a.i == p))
    # Upper flow bound at T nodes
    @constraint(pooling_problem, [t in T], sum(x[a] for a in A if a.j == t) <= N[t].bub)
    # Sulfur balance at P nodes
    @constraint(pooling_problem, [p in P], sum(q[a.i]*x[a] for a in A if a.j == p) == q[p]*sum(x[a] for a in A if a.i == p))
    # Sulfur balances at T nodes
    @constraint(pooling_problem, [t in T], sum(q[a.i]*x[a] for a in A if a.j == t) == q[t]*sum(x[a] for a in A if a.j == t))
    # Sulfur upper bounds at T nodes
    @constraint(pooling_problem, [t in T], q[t] <= N[t].qub)
    # Sulfur values at S nodes (here N[s].qlb == N[s].qub so either can be used)
    @constraint(pooling_problem, [s in S], q[s] == N[s].qub)

    # set start values for pool nodes
    for i = 1:length(P)
        set_start_value(pooling_problem[:q][P[i]], pool_start_value[i])
    end

    # Print pooling_problem. NOTE: Print pooling_problem at any point to see how it looks
    return pooling_problem, bub[findall(x -> x != Inf, bub)]
end 

