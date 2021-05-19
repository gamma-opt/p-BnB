"""
    penalty_parameter_update(initial_parameters::MIP_initial_parameters, generated_parameters::MIP_generated_parameters, x::Array{Float64})
returns the updated value for each component of the penalty parameter vector
"""

function penalty_parameter_update(initial_parameters::MIP_initial_parameters, generated_parameters::MIP_generated_parameters, x::Array{Float64})

    # calculating the exepected value
    x_avg = round.(sum(initial_parameters.scen_prob[j] .* x[:,j] for j = 1:initial_parameters.num_scen), digits = 8)

    #calculating x max value among scenarios
    x_max = findmax(x, dims = 2)[1]

    # calculating x min value among scenarios
    x_min = findmin(x, dims = 2)[1]

    # defining the array for penalty values
    ρ = zeros(initial_parameters.num_first_stage_var)

    for i_int in generated_parameters.x_int_indexes
        ρ[i_int] = generated_parameters.objective_c[i_int] / (x_max[i_int] - x_min[i_int] + 1)
    end

    for i_cont in generated_parameters.x_cont_indexes
        ρ[i_cont] = generated_parameters.objective_c[i_cont] / max( sum( initial_parameters.scen_prob[j] *  abs(x[i_cont, s] - x_avg[i_cont]) for s = 1:initial_parameters.num_scen ) , 1)
    end

     return ρ

end
