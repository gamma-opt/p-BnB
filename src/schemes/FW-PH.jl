function FW_PH( bnb_node::node, V_0::Array{AbstractArray{Vector{Array{Float64}},1}}, x_0::Array{Array{Float64}}, w_0::Array{Array{Float64}} )

    # simplifying the notations
    initial_parameters = bnb_node.initial_parameters
    generated_parameters = bnb_node.generated_parameters

    α = initial_parameters.PH_SDM_parameters.α
    PH_ϵ = initial_parameters.PH_SDM_parameters.PH_tol
    k_max = initial_parameters.PH_SDM_parameters.PH_max_iter
    SDM_ϵ = initial_parameters.PH_SDM_parameters.SDM_tolerance
    t_max = initial_parameters.PH_SDM_parameters.SDM_max_iter

    # calculating starting value for parameter z
    z_0 = sum(initial_parameters.scen_prob[s] .* x_0[s]  for s = 1:initial_parameters.num_scen)
    #@show z_0
    #@show x_0[1]

    w_k =  Array{Array{Float64}}(undef, k_max, initial_parameters.num_scen)
    [ w_k[1,s] = w_0[s] .+ initial_parameters.al_penalty_parameter .* (x_0[s] .- z_0) for s = 1:initial_parameters.num_scen ]

    z_k = Array{Array{Float64}}(undef, k_max)

    # variable to store the dual value
    ϕ_k = Array{Float64}(undef, k_max)

    # defining the variables for saving the results of SDM
    x_k = Array{Array{Float64}}(undef, k_max, initial_parameters.num_scen)
    y_k = Array{Array{Float64}}(undef, k_max, initial_parameters.num_scen)
    w_RNMDT_k = Array{Array{Float64}}(undef, k_max, initial_parameters.num_scen)
    z_FR_k = Array{Array{Float64}}(undef, k_max, initial_parameters.num_scen)
    V_k = Array{AbstractArray{Vector{Array{Float64}},1}}(undef,k_max, initial_parameters.num_scen)
    #V_k = Array{Any}(undef, k_max, initial_parameters.num_scen)
    dual_value_k = Array{Float64}(undef, k_max, initial_parameters.num_scen)

    for k = 1:k_max
    #for k = 1:1

        #print("\n\n\nFW-PH ITERATION $k\n\n\n")
        @sync Base.Threads.@threads for s = 1:initial_parameters.num_scen
        #@show s
        #for s = 1:initial_parameters.num_scen
        #print("\n\n\nSCENARIO $s\n\n\n")

            if k == 1

                x_wave_s = (1 - α) .* z_0 + α .* x_0[s]
                x_k[k,s], y_k[k,s], w_RNMDT_k[k,s], z_FR_k[k,s], V_k[k,s], dual_value_k[k,s] = SDM(s, bnb_node, V_0[s], x_wave_s, V_0[s][1][2], V_0[s][1][3], V_0[s][1][4], w_k[k,s], z_0, t_max, SDM_ϵ)
                #@show x_k[k,s]
            else

                x_wave_s = (1 - α) .* z_k[k-1] .+ α .* x_k[k-1,s]
                x_k[k,s], y_k[k,s], w_RNMDT_k[k,s], z_FR_k[k,s], V_k[k,s], dual_value_k[k,s] = SDM(s, bnb_node, V_k[k-1,s], x_wave_s, y_k[k-1,s], w_RNMDT_k[k-1,s], z_FR_k[k-1,s], w_k[k,s], z_k[k-1], t_max, SDM_ϵ)


            end
        end

        ϕ_k[k] = sum(initial_parameters.scen_prob[s] .* dual_value_k[k,s]  for s = 1:initial_parameters.num_scen)

        z_k[k] = sum(initial_parameters.scen_prob[s] .* x_k[k,s]  for s = 1:initial_parameters.num_scen)

        #@show ϕ_k[k]
        #@show z_k[k]
        #@show x_k[k,:]

        if sqrt(sum( initial_parameters.scen_prob[s] .* norm(x_k[k,s] .- (k == 1 ? z_0 : z_k[k-1]) )^2 for s = 1 : initial_parameters.num_scen)) < PH_ϵ
            return(x_k[k, :], y_k[k, :], w_RNMDT_k[k, :], z_FR_k[k,:], z_k[k], w_k[k, :], ϕ_k[k])
        end

        #@show k

        if k < k_max
            [ w_k[k+1, s] = w_k[k, s] .+ initial_parameters.al_penalty_parameter .* (x_k[k,s] .- z_k[k]) for s = 1: initial_parameters.num_scen]
        end

    end



    #return(x_k[1,:], y_k[1,:], w_RNMDT_k[1,:])
    return (x_k[end, :], y_k[end, :], w_RNMDT_k[end, :], z_FR_k[end,:], z_k[end], w_k[end, :], ϕ_k[end])

end
