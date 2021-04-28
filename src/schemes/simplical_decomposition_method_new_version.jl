"""
    SDM(V_t::Array{Array{Array{Float64}}}, )
returns an mixed-integer optimisation problem
a JuMP::model with the parameters token from "intial_parameters"

"""
function SDM_new(scenario::Int, bnb_node::node, V_0::Array{Array{Float64}}, x_0::Array{Float64}, w_s::Array{Float64}, z_SDM::Array{Float64}, t_max::Int, τ::Float64)

    # variable for storing the dual value
    dual_value_s = 0

    # formulating initial finite set of points
    V_t = Array{Any}(undef, t_max+1)
    V_t[1] = V_0

    # simplifying the notations
    initial_parameters = bnb_node.initial_parameters
    generated_parameters = bnb_node.generated_parameters

    # variable to store updated lagrangian multipliers values
    w_t = Array{Array{Float64}}(undef, t_max)

    # variable to store updated first-stage variables values
    x_t = Array{Array{Float64}}(undef, t_max)

    # variable to store updated second-stage variables values
    y_t = Array{Array{Float64}}(undef, t_max)

    # variable to store updated RNMDT-auxiliary variables w values
    w_RNMDT_t = Array{Array{Float64, 2}}(undef, t_max)

    # bound gap meassure
    Γ_t = Array{Float64}(undef, 1, t_max)

    # auxiliary augmented lagrangian approximation problem
    al_approximation = copy(bnb_node.dual_subproblems[scenario])
    set_optimizer(al_approximation, optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol, "Method" => initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag, "Threads" => initial_parameters.gurobi_parameters.Threads, "NumericFocus" => initial_parameters.gurobi_parameters.NumericFocus))

    for t = 1:t_max

        w_t[t] = w_s .+ initial_parameters.al_penalty_parameter.* ((t==1 ? x_0 : x_t[t-1]) .- z_SDM)

        @objective( al_approximation, Min,
            -
            ( sum(generated_parameters.objective_Qs[scenario][i, j] * al_approximation[:w_RNMDT][i, j]
                for i = 1 : initial_parameters.num_second_stage_var,
                    j = 1 : initial_parameters.num_second_stage_var)
                + sum( al_approximation[:x][i] * generated_parameters.objective_c[i]  for i = 1:initial_parameters.num_first_stage_var)
                + sum( al_approximation[:y][j] * generated_parameters.objective_fs[scenario][j]  for j = 1:initial_parameters.num_second_stage_var)

                - sum( w_t[t] .* al_approximation[:x] )
                - initial_parameters.μ * sum(al_approximation[:z][r] for r  = 1:initial_parameters.num_const )

            )

        )

        status = optimize!(al_approximation)

        # variable to store updated approximated first-stage variables values
        x_hat = value.(al_approximation[:x])

        # variable to store updated approximated second-stage variables values
        y_hat =  value.(al_approximation[:y])

        # variable to store updated approximated RNMDT-auxiliary variables w values
        w_RNMDT_hat = value.(al_approximation[:w_RNMDT])

        # if we are at the very first iteration we use the starting values
        if t == 1

            # updating dual value
            dual_value_s = objective_value(al_approximation)

            # calculating the value of the bound gap at iteration t == 1
            Γ_t_value  =  sum(generated_parameters.objective_Qs[scenario][i, j] * (w_RNMDT_hat[i, j] - V_0[3][i,j])
                for i = 1 : initial_parameters.num_second_stage_var,
                    j = 1 : initial_parameters.num_second_stage_var)
                + sum( (x_hat[i] - x_0[i]) * generated_parameters.objective_c[i]  for i = 1:initial_parameters.num_first_stage_var)
                + sum( (y_hat[j] - V_0[2][j]) * generated_parameters.objective_fs[scenario][j]  for j = 1:initial_parameters.num_second_stage_var)
                - sum( w_t[t] .* (x_hat .- x_0) )
                - initial_parameters.μ * sum(value.(al_approximation[:z])[r] for r  = 1:initial_parameters.num_const )
        else
            # calculating the value of the bound gap at iteration t
            Γ_t_value  =  sum(generated_parameters.objective_Qs[scenario][i, j] * (w_RNMDT_hat[i, j] - w_RNMDT_t[t-1][i,j])
                for i = 1 : initial_parameters.num_second_stage_var,
                    j = 1 : initial_parameters.num_second_stage_var)
                + sum( (x_hat[i] - x_t[t-1][i]) * generated_parameters.objective_c[i]  for i = 1:initial_parameters.num_first_stage_var)
                + sum( (y_hat[j] - y_t[t-1][j]) * generated_parameters.objective_fs[scenario][j]  for j = 1:initial_parameters.num_second_stage_var)
                - sum( w_t[t] .* (x_hat .- x_t[t-1]) )
                - initial_parameters.μ * sum(value.(al_approximation[:z])[r] for r  = 1:initial_parameters.num_const )

        end

        # adding new values for bound gap and finite set of points at iteration t
        Γ_t[t] = Γ_t_value
        V_t[t+1] = [x_hat, y_hat, w_RNMDT_hat] # since we also had 0 element

        # formulating new problem representing augmented lagrangian
        al_SDM = copy(bnb_node.dual_subproblems[scenario])

        @suppress set_optimizer(al_SDM, optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol, "Method" => initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag,  "Threads" => initial_parameters.gurobi_parameters.Threads, "NumericFocus" => initial_parameters.gurobi_parameters.NumericFocus, "Presolve" => 0))

        # fixing the value of the variable al_z
        fix.(al_SDM[:al_z], z_SDM)
        JuMP.unset_integer.(al_SDM[:al_z])

        # defining current length of V_t
        cl_V_t = t+1 # since we have 0 elemnt as well - starting values

        # defining new variables and constraints
        @variable(al_SDM, a[1:cl_V_t]>=0)
        @constraint(al_SDM, sum(a)==1)

        @constraint(al_SDM, al_SDM[:x] .== sum(a[i] .* V_t[i][1] for i = 1:cl_V_t ))
        @constraint(al_SDM, al_SDM[:y] .== sum(a[i] .* V_t[i][2] for i = 1:cl_V_t ))
        @constraint(al_SDM, al_SDM[:w_RNMDT] .== sum(a[i] .* V_t[i][3] for i = 1:cl_V_t ))

        #print(al_SDM)
        # optimising problem representing augmented lagrangian
        @suppress optimize!(al_SDM)

        # updating the values of the primal variables
        x_t[t] = value.(al_SDM[:x])
        y_t[t] = value.(al_SDM[:y])
        w_RNMDT_t[t] = value.(al_SDM[:w_RNMDT])

        #if bound gap is smaller than predefiend tolerance
        if Γ_t_value <= τ
            return(x_t[t], y_t[t], w_RNMDT_t[t], V_t[t+1], dual_value_s, Γ_t[t])
        end

    end

    return(x_t[end], y_t[end], w_RNMDT_t[end], V_t[end], dual_value_s, Γ_t[end] )

end
