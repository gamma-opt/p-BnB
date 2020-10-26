"""
    dp_based_RNMDT_input
Stores the parameters for dynamic precision based RNMDT method. Has the following fields:
* `N1_percentage::Int64`:                   Perecnatuge that is used to calculate N1 as some percent of the total number
                                            of second stage decision variables in the primal problem
* `N2::Int64`:                              iterations related parameter
* `tolerance::Float64`:                     tolerance used for the algorithm (while (UB_optimal - LB_optimal)/LB_optimal > tolerance)
* `time_limit::Float64`:                    time limit for the algorithm
* `max_number_of_iterations::Int64`:        iterations limit for the alogrothm

"""
mutable struct dp_based_RNMDT_input
    N1_percentage::Int64
    N2::Int64
    tolerance::Float64
    time_limit::Float64
    max_number_of_iterations::Int64
end

"""
    dp_based_RNMDT_ouput
    bm_output
Stores the output of the dynamic precision-based RNMDT method. Has the following fields:
* `dual_objective_value::Vector{Float64}`:        Contains the values of the dual objective function for all the iterations
* `variables_values::Array{Any}`:                 Contains the values of the first-stage decision variables at the final iteration as the first element,
                                                           the values of the second-stage decision variables at the final iteration as the second element,

"""
mutable struct dp_based_RNMDT_output
    dual_objective_value::Array{Float64}
    variables_values::Array{Array{Float64}}
end
