function FW_PH_V_0_initialisation( bnb_node::node, w_s::Array{Array{Float64}})

    # simplifying the notation
    initial_parameters = bnb_node.initial_parameters
    generated_parameters = bnb_node.generated_parameters
    dual_subproblems = Array{JuMP.Model}(undef, initial_parameters.num_scen)

    x_0 = Array{Array{Float64}}(undef, initial_parameters.num_scen)
    y_0 = Array{Array{Float64}}(undef, initial_parameters.num_scen)
    w_RNMDT_0 = Array{Array{Float64}}(undef, initial_parameters.num_scen)
    z_0 = Array{Array{Float64}}(undef, initial_parameters.num_scen)

    y_1 = Array{Array{Float64}}(undef, initial_parameters.num_scen)
    w_RNMDT_1 = Array{Array{Float64}}(undef, initial_parameters.num_scen)
    z_1 = Array{Array{Float64}}(undef, initial_parameters.num_scen)

    V_0 = Array{AbstractArray{Vector{Array{Float64}},1}}(undef, initial_parameters.num_scen)

    for s = 1:initial_parameters.num_scen
        dual_subproblems[s] = copy(bnb_node.dual_subproblems[s])
        @suppress set_optimizer(dual_subproblems[s], optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol, "Method" => initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag,        "Threads" => initial_parameters.gurobi_parameters.Threads, "NumericFocus" => initial_parameters.gurobi_parameters.NumericFocus, "Presolve" => 0))

        @objective(dual_subproblems[s], Min,
            -
            ( sum(generated_parameters.objective_Qs[s][i, j] * dual_subproblems[s][:w_RNMDT][i, j]
                for i = 1 : initial_parameters.num_second_stage_var,
                    j = 1 : initial_parameters.num_second_stage_var)
                + sum( dual_subproblems[s][:x][i] * generated_parameters.objective_c[i]  for i = 1:initial_parameters.num_first_stage_var)
                + sum( dual_subproblems[s][:y][j] * generated_parameters.objective_fs[s][j]  for j = 1:initial_parameters.num_second_stage_var)
            )
           + sum( w_s[s] .* (dual_subproblems[s][:x]) )

           + initial_parameters.μ * sum(dual_subproblems[s][:z][r] for r  = 1:initial_parameters.num_const )
        )

        @suppress optimize!(dual_subproblems[s])

        x_0[s] = value.(dual_subproblems[s][:x])
        y_0[s] = value.(dual_subproblems[s][:y])
        w_RNMDT_0[s] = value.(dual_subproblems[s][:w_RNMDT])
        z_0[s] = value.(dual_subproblems[s][:z])

        #@show z_0[s]
        #@show x_0[s]

        if s == 1
            V_0[s] = [[ x_0[s], y_0[s], w_RNMDT_0[s], z_0[s] ]]
        else
            fix.(dual_subproblems[s][:x], x_0[1])
            @objective(dual_subproblems[s], Min,
                -
                ( sum(generated_parameters.objective_Qs[s][i, j] * dual_subproblems[s][:w_RNMDT][i, j]
                    for i = 1 : initial_parameters.num_second_stage_var,
                        j = 1 : initial_parameters.num_second_stage_var)
                    + sum( dual_subproblems[s][:x][i] * generated_parameters.objective_c[i]  for i = 1:initial_parameters.num_first_stage_var)
                    + sum( dual_subproblems[s][:y][j] * generated_parameters.objective_fs[s][j]  for j = 1:initial_parameters.num_second_stage_var)
                )

               + initial_parameters.μ * sum(dual_subproblems[s][:z][r] for r  = 1:initial_parameters.num_const )
               )

             #print(dual_subproblems[s])
             @suppress optimize!(dual_subproblems[s])
             y_1[s] = value.(dual_subproblems[s][:y])
             w_RNMDT_1[s] = value.(dual_subproblems[s][:w_RNMDT])
             z_1[s] = value.(dual_subproblems[s][:z])
             V_0[s] = [[ x_0[s], y_0[s], w_RNMDT_0[s], z_0[s] ], [ x_0[1], y_1[s], w_RNMDT_1[s], z_1[s]]]
        end
    end


    return V_0, x_0

end
