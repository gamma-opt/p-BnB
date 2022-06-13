function bounds_tightening(initial_parameters::MIP_initial_parameters,generated_parameters::MIP_generated_parameters, McCromick_relaxation::JuMP.Model, f0::Float64)
    new_y_boundaries = Array{Any}(undef, initial_parameters.num_scen)
    auxiliary_array = Array{Float64}(undef,initial_parameters.num_second_stage_var, 2)
    #new_y_boundaries = Array{Float64}(undef, initial_parameters.num_second_stage_var, initial_parameters.num_scen, 2)

    bounds_tightening_problem = copy(McCromick_relaxation)
    @suppress set_optimizer(bounds_tightening_problem, optimizer_with_attributes(Gurobi.Optimizer,  "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol,
    "Method" => initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag, "Threads" => initial_parameters.gurobi_parameters.Threads))
    JuMP.unset_integer.(bounds_tightening_problem[:x][:, :])

    @constraint(bounds_tightening_problem,
        sum( initial_parameters.scen_prob[s] *
            (
                sum( bounds_tightening_problem[:y][i, s] * generated_parameters.objective_Qs[s][i, j] * bounds_tightening_problem[:y][j, s] for i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var)
                + sum( bounds_tightening_problem[:x][i, s] * generated_parameters.objective_c[i] for i = 1 : initial_parameters.num_first_stage_var)
                + sum( bounds_tightening_problem[:y][i, s] * generated_parameters.objective_fs[s][i] for i = 1 : initial_parameters.num_second_stage_var)
                
            )
            + initial_parameters.μ * sum( bounds_tightening_problem[:z][s,r] for r  = 1:initial_parameters.num_const )
    for s in 1:initial_parameters.num_scen)
        <= f0
    )

    #print(bounds_tightening_problem)
    for s = 1:initial_parameters.num_scen
        for i = 1:initial_parameters.num_second_stage_var

            @objective(bounds_tightening_problem, Min, bounds_tightening_problem[:y][i,s] )
            optimize!(bounds_tightening_problem)
            auxiliary_array[i,1] = objective_value(bounds_tightening_problem)

            @objective(bounds_tightening_problem, Max, bounds_tightening_problem[:y][i,s] )
            optimize!(bounds_tightening_problem)
            auxiliary_array[i,2] = objective_value(bounds_tightening_problem)
            @show auxiliary_array
        end

        new_y_boundaries[s] = copy(auxiliary_array)
        @show new_y_boundaries

    end

    return new_y_boundaries

end
