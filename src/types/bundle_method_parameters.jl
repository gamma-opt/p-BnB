"""
    bm_input
Stores the input parameters of the bundle method. Has the following fields:

* `parallelisation_is_used::Bool`:                                  Is the alogrithm parallelised scenario wise (true = yes, false = no)

* `max_number_of_iterations::Int64`:                                Maximum number of the iterations
* `number_of_iteration_for_checking::Int64`:                        Parameter indicating for how many iterations we check the stopping criteria
* `eps_stop::Float64`:                                              stoping criteria tolerance

* `m::Int64`:                                                       Acceptance parameter for the serious step condition
* `d::Int64`:                                                       Stepsize parameter for the bundle method
* `initial_centre_of_gravity_min::Float64`:                         Minimum value for the initial centre of mass
* `initial_centre_of_gravity_max::Float64`:                         Maximum value for the initial centre of mass
"""


mutable struct bm_input

    # general parameters
    parallelisation_is_used::Bool

    # stoping criteria parameters
    max_number_of_iterations::Int64
    number_of_iteration_for_checking::Int64
    eps_stop::Float64

    # algorithm performace parameters
    m::Float64
    d::Float64
    initial_centre_of_gravity_min::Float64
    initial_centre_of_gravity_max::Float64

end


"""
    bm_output
Stores the output of the bundle method. Has the following fields:
* `dual_objective_value::Vector{Float64}`:        Contains the values of the dual objective function for all the iterations
* `variables_values::Array{Any}`:                 Contains the values of the first-stage decision variables at the final iteration as the first element,
                                                           the values of the second-stage decision variables at the final iteration as the second element,
                                                            * if bundle method is applied to RNMDT-based relaxation
                                                                the values of the auxiliary variable w at the final iteration as the third element
                                                                the values of centre of mass at the final iteration as the fourth element
"""

struct bm_output
    dual_objective_value::Array{Float64}
    variables_values::Array{Array{Float64}}
end
