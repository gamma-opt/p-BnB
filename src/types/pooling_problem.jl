## This file contains all the necessary types required for pooling problems 

##---------------Types to define the pooling problem-----------------------
## Struct for arcs 
struct Arc
    i::String      # Head node of arc (i,j) ∈ A
    j::String      # Tail node of arc (i,j) ∈ A
end

## Struct for nodes
struct Node
    cost::Float64  # Node costs for s∈S and profits for t∈T
    qlb::Float64   # Lower property value bound
    qub::Float64   # Upper property value bound
    bub::Float64   # Upper flow bound
end

##---------------Types to extract the infomation about the constraints----
##----and objective functions coefficents from the pooling problem--------


# structure to remember each variable corresponds to each index
struct pooling_variable
    name::String
    index::Int64
end

# strucuture to keep the coefficeints for quadratic affine constraint
struct quad_const
    quadratic_coefficients_matrix::Array{Float64}
    affine_coefficients_vector::Array{Float64}
    sign::String
    constant::Float64
end

# strucuture to keep coefficents for affine constraint
struct af_const
    linear_coefficients_matrix::Array{Float64}
    sign::String
    RHS::Array{Float64}
end 

# strucuture to keep one-sided of the box constraints, i.e., x<=c or x>= c or x==c 
struct box_const
    variable_index::Int64
    sign::String
    RHS::Float64
end 

##---------------Types needed to define modified stochastic-----------
##---------------two-stage pooling problem----------------------------

struct stochastic_pooling_problem_parameters
    pool_max_capacity::Float64
    flow_max_capacity::Float64

    demand_min_value::Float64
    demand_max_value::Float64

    source_node_min_value::Float64
    source_node_max_value::Float64
end

##---------------Structure needed to contain all parameters-----------
##---------------related to pooling problem---------------------------

struct pooling_problem_parameters
    primal_pool_problem_link::String
    print_out_variables_correlation_map::Bool
    start_value_pooling_nodes::Array{Float64}
    flow_cost::Float64
    pool_cost::Float64
    stoc_pp::stochastic_pooling_problem_parameters
end