 function al_dual_value_calculation(initial_parameters::MIP_initial_parameters, generated_parameters:: MIP_generated_parameters, x::Array{Array{Float64}}, y::Array{Array{Float64}}, w_RNMDT::Array{Array{Float64}}, al_z::Array{Float64}, vector_of_lambda_lagrangian::Array{Array{Float64}})
     al_dual_value = Array{Float64}(undef, initial_parameters.num_scen)
     for s = 1:initial_parameters.num_scen
         al_dual_value[s] =  - initial_parameters.scen_prob[s] *
                     ( sum(generated_parameters.objective_Qs[s][i, j] * w_RNMDT[s][i, j]
                         for i = 1 : initial_parameters.num_second_stage_var,
                             j = 1 : initial_parameters.num_second_stage_var)
                     + sum( x[s][i] * generated_parameters.objective_c[i]  for i = 1:initial_parameters.num_first_stage_var)
                     + sum( y[s][j] * generated_parameters.objective_fs[s][j]  for j = 1:initial_parameters.num_second_stage_var)

                     + (1/initial_parameters.scen_prob[s]) * sum( vector_of_lambda_lagrangian[s] .* (x[s] .- al_z) )
                     + initial_parameters.al_penalty_parameter/2 * sum( (x[s] .- al_z) .* (x[s] .- al_z) )

                     #- initial_parameters.μ * sum(z[r] for r  = 1:initial_parameters.num_const )

                     )
     end
     return(sum(al_dual_value))

 end
 function FW_PH_bundle( bnb_node::node, V_0::Array{Any}, x_0::Array{Float64,2}, y_0::Array{Float64,2}, w_RNMDT_0::Array{Float64,3}, w_0::Array{Array{Float64}} )

    γ = 0.7
    # simplifying the notations
    initial_parameters = bnb_node.initial_parameters
    generated_parameters = bnb_node.generated_parameters

    α = initial_parameters.PH_SDM_parameters.α
    PH_ϵ = initial_parameters.PH_SDM_parameters.PH_tol
    k_max = initial_parameters.PH_SDM_parameters.PH_max_iter
    SDM_ϵ = initial_parameters.PH_SDM_parameters.SDM_tolerance
    t_max = initial_parameters.PH_SDM_parameters.SDM_max_iter

    # calculating starting value for parameter z
    z_0 = sum(initial_parameters.scen_prob[s] .* x_0[:,s]  for s = 1:initial_parameters.num_scen)

    w_k =  Array{Array{Float64}}(undef, k_max, initial_parameters.num_scen)
    #[ w_k[1,s] = w_0[s] .+ initial_parameters.al_penalty_parameter .* (x_0[s] .- z_0) for s = 1:initial_parameters.num_scen ]

    z_k = Array{Array{Float64}}(undef, k_max)

    # variable to store the dual value
    ϕ_k = Array{Float64}(undef, k_max)

    # defining the variables for saving the results of SDM
    x_k = Array{Array{Float64}}(undef, k_max, initial_parameters.num_scen)
    y_k = Array{Array{Float64}}(undef, k_max, initial_parameters.num_scen)
    w_RNMDT_k = Array{Array{Float64}}(undef, k_max, initial_parameters.num_scen)
    V_k = Array{Any}(undef, k_max, initial_parameters.num_scen)
    dual_value_k = Array{Float64}(undef, k_max, initial_parameters.num_scen)
    Γ_value_k = Array{Float64}(undef, k_max, initial_parameters.num_scen)
    γ_value_k = Array{Float64}(undef, k_max)

    ## zero step
    Γ_value_zero_iter = Array{Float64}(undef, initial_parameters.num_scen)
    dual_value_zero_iter = Array{Float64}(undef, initial_parameters.num_scen)
    x_zero_iter = Array{Array{Float64}}(undef, initial_parameters.num_scen)
    y_zero_iter = Array{Array{Float64}}(undef, initial_parameters.num_scen)
    w_RNMDT_zero_iter = Array{Array{Float64}}(undef, k_max, initial_parameters.num_scen)
    for s = 1:initial_parameters.num_scen
        x_wave_s = (1 - α) .* z_0 + α .* x_0[:,s]
        x_zero_iter[s], y_zero_iter[s], w_RNMDT_zero_iter[s], V_0[s], dual_value_zero_iter[s], Γ_value_zero_iter[s] = SDM_new(s, bnb_node, V_0[s], x_wave_s, w_0[s], z_0, t_max, SDM_ϵ)
    end

    z_0 = sum(initial_parameters.scen_prob[s] .* x_zero_iter[:,s]  for s = 1:initial_parameters.num_scen)

    ϕ_waveline = al_dual_value_calculation(initial_parameters, generated_parameters, x_zero_iter, y_zero_iter, w_RNMDT_zero_iter, z_0, w_0) + initial_parameters.al_penalty_parameter/2 * sum( sum( (x_zero_iter[s] .- z_0) .* (x_zero_iter[s] .- z_0) ) for s = 1:initial_parameters.num_scen) - sum(Γ_value_zero_iter)

    ϕ_underline_k = Array{Float64}(undef, k_max)
    ϕ_underline_0 = ϕ_waveline


    # initialisation for the very first iteraion
    ϕ_underline_k[1] = ϕ_underline_0
    [ w_k[1,s] = w_0[s] .+ initial_parameters.al_penalty_parameter .* (x_zero_iter[s] .- z_0) for s = 1:initial_parameters.num_scen ]

    for k = 1:k_max

        # initialisation
        if k > 1
            ϕ_underline_k[k] = ϕ_underline_k[k-1]
            w_k[k,:] = w_k[k-1,:]
        end

        for s = 1:initial_parameters.num_scen
        #@show s

            if k == 1

                x_wave_s = (1 - α) .* z_0 + α .* x_zero_iter[s]
                x_k[k,s], y_k[k,s], w_RNMDT_k[k,s], V_k[k,s], dual_value_k[k,s], Γ_value_k[k,s] = SDM_new(s, bnb_node, V_0[s], x_wave_s, w_k[k,s], z_0, t_max, SDM_ϵ)

            else

                x_wave_s = (1 - α) .* z_k[k-1] .+ α .* x_k[k-1,s]
                x_k[k,s], y_k[k,s], w_RNMDT_k[k,s], V_k[k,s], dual_value_k[k,s], Γ_value_k[k,s] = SDM_new(s, bnb_node, V_k[k-1,s], x_wave_s, w_k[k,s], z_k[k-1], t_max, SDM_ϵ)


            end
        end

        ϕ_k[k] = sum(initial_parameters.scen_prob[s] .* dual_value_k[k,s]  for s = 1:initial_parameters.num_scen)
        z_k[k] = sum(initial_parameters.scen_prob[s] .* x_k[k,s]  for s = 1:initial_parameters.num_scen)


        if sqrt(sum( initial_parameters.scen_prob[s] .* norm(x_k[k,s] .- z_k[k])^2 for s = 1 : initial_parameters.num_scen)) < PH_ϵ
            return(x_k[k, :], y_k[k, :], w_RNMDT_k[k, :], z_k[k], w_k[k, :], ϕ_k[k])
        end

        ϕ_waveline = al_dual_value_calculation(initial_parameters, generated_parameters, x_k[k,:], y_k[k,:], w_RNMDT_k[k,:], z_k[k], w_k[k,:])
                     + initial_parameters.al_penalty_parameter/2 * sum( sum( (x_k[k,s] .- z_0) .* (x_k[k,s] .- z_0) ) for s = 1:initial_parameters.num_scen) - sum( Γ_value_k[k,:] )

        γ_value_k[k] = (ϕ_waveline - ϕ_underline_k[k]) / ( al_dual_value_calculation(initial_parameters, generated_parameters, x_k[k,:], y_k[k,:], w_RNMDT_k[k,:], z_k[k], w_k[k,:]) + initial_parameters.al_penalty_parameter/2 * sum( sum( (x_k[k,s] .- z_0) .* (x_k[k,s] .- z_0) ) for s = 1:initial_parameters.num_scen) - ϕ_underline_k[k])

        if γ_value_k[k] >= γ
            [ w_k[k, s] = w_k[k, s] + initial_parameters.al_penalty_parameter .* (x_k[k,s] .- z_k[k]) for s = 1: initial_parameters.num_scen]
            ϕ_underline_k[k] = ϕ_waveline
        end

    end



    #return(x_k[1,:], y_k[1,:], w_RNMDT_k[1,:])
    return (x_k[end, :], y_k[end, :], w_RNMDT_k[end, :], z_k[end], w_k[end, :], ϕ_k[end])

end
