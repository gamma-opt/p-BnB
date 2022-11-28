"""
    MIP_generation(intial_parameters::MIP_initial_parameters)
returns an mixed-integer optimisation problem
a JuMP::model with the parameters token from "intial_parameters"

"""

function MIP_generation(initial_parameters::MIP_initial_parameters, generated_parameters::MIP_generated_parameters)

    # depending on whether the variables are fixed or not we use Ipopt and Gurobi as a solver respectively
    original_problem = Model( optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol, "Method" => initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag,  "Threads" => initial_parameters.gurobi_parameters.Threads, "TimeLimit" => initial_parameters.solver_time_limit ) )

    # first stage decision variables
    @variable(original_problem, x[ 1 : initial_parameters.num_first_stage_var, 1 : initial_parameters.num_scen ], Int )
    JuMP.unset_integer.(original_problem[:x][generated_parameters.x_cont_indexes, :])

    # second stage decision variables
    @variable(original_problem, y[ 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_scen ], Int)
    [JuMP.unset_integer.(original_problem[:y][generated_parameters.y_cont_indexes[:,j], j]) for j = 1:initial_parameters.num_scen]

    # slack variables
    @variable(original_problem, z[ 1 : initial_parameters.num_scen, 1 : initial_parameters.num_const ] >=0 )

    # quadratic objective
    @objective(original_problem, Min,
        #map_coefficients_inplace!(a -> round(a, digits=1),
        -sum( initial_parameters.scen_prob[s] *
            (
                sum( y[i, s] * generated_parameters.objective_Qs[s][i, j] * y[j, s] for i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var)
                + sum( x[i, s] * generated_parameters.objective_c[i] for i = 1 : initial_parameters.num_first_stage_var)
                + sum( y[i, s] * generated_parameters.objective_fs[s][i] for i = 1 : initial_parameters.num_second_stage_var)
                - initial_parameters.μ * sum(z[s,r] for r  = 1:initial_parameters.num_const )
            )
        for s in 1:initial_parameters.num_scen)
        )
        #)

    # quadratic constraints
    @constraint(original_problem, [s = 1 : initial_parameters.num_scen, i = 1 : initial_parameters.num_const ],
        sum( y[j, s] * generated_parameters.constraint_Qs[s][i][j,k] * y[k, s] for j = 1 : initial_parameters.num_second_stage_var, k = 1: initial_parameters.num_second_stage_var)
        + sum( x[j, s] * generated_parameters.constraint_fs[s][i][1, j] for j = 1 : initial_parameters.num_first_stage_var)
        + sum( y[j, s] * generated_parameters.constraint_fs[s][i][2, j] for j = 1:  initial_parameters.num_second_stage_var)
        + generated_parameters.constraint_fs[s][i][3, 1] + z[s, i]<= 0 )

    # box constraints for first stage variables
    @constraint(original_problem, [ s = 1 : initial_parameters.num_scen ],
        generated_parameters.x_boundaries[:, 1] .<= x[:, s] .<= generated_parameters.x_boundaries[:, 2])

    # box constraints for second stage variables
    @constraint(original_problem, [s = 1 : initial_parameters.num_scen],
        generated_parameters.y_boundaries[s][:, 1] .<= y[:, s] .<= generated_parameters.y_boundaries[s][:, 2])

    # non-anticipativity conditions
    @constraint( original_problem, [s in 2 : initial_parameters.num_scen, i = 1 : initial_parameters.num_first_stage_var],
        x[i, s] - x[i, 1] == 0 )

    return original_problem

end

"""
function f_lambda_lagrangian(lambda_lagrangian::Array{Any}, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1:end]) : - lambda_lagrangian[dec_index-1])

"""

# auxiliary function for Lagrangian multipliers update
f_lambda_lagrangian(lambda_lagrangian::Array{Array{Float64}}, dec_index::Int) = (dec_index == 1 ? sum(lambda_lagrangian[1:end]) : - lambda_lagrangian[dec_index-1])

function primal_problem_based_lagrangian_relaxation_generation(initial_parameters::MIP_initial_parameters, generated_parameters::MIP_generated_parameters)

    # generate the starting values for the lagrangian multipliers
    vector_of_lambda_lagrangian = Array{Array{Float64}}(undef, initial_parameters.num_scen - 1)
    [ vector_of_lambda_lagrangian[i] = zeros(1,initial_parameters.num_first_stage_var)
            for i = 1 : initial_parameters.num_scen - 1 ]

    # creating the array of subproblems
    vector_of_subproblems = Array{Any}(undef, initial_parameters.num_scen)

    # formulating the subproblems based on the scenarios
    for s = 1 : initial_parameters.num_scen

        vector_of_subproblems[s] = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol, "Method" => initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag,  "Threads" => initial_parameters.gurobi_parameters.Threads, "NumericFocus" => initial_parameters.gurobi_parameters.NumericFocus))

        # first stage decision variables
        @variable( vector_of_subproblems[s], x[1 : initial_parameters.num_first_stage_var], Int )
        JuMP.unset_integer.(vector_of_subproblems[s][:x][generated_parameters.x_cont_indexes])

        # second stage decision variables
        @variable( vector_of_subproblems[s], y[1 : initial_parameters.num_second_stage_var], Int )
        JuMP.unset_integer.(vector_of_subproblems[s][:y][generated_parameters.y_cont_indexes[:,s]])

        # slack variables
        @variable(vector_of_subproblems[s], z[ 1 : initial_parameters.num_const ] >=0 )

        # quadratic objective
        @objective( vector_of_subproblems[s], Min,
          - initial_parameters.scen_prob[s] *
            (
            sum( y[i] * initial_parameters.scen_prob[s] * generated_parameters.objective_Qs[s][i, j] * y[j] for i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var)
            + sum( x[i] * initial_parameters.scen_prob[s] * generated_parameters.objective_c[i]  for i = 1:initial_parameters.num_first_stage_var)
            + sum( y[j] * initial_parameters.scen_prob[s] * generated_parameters.objective_fs[s][j]  for j = 1:initial_parameters.num_second_stage_var)
            - initial_parameters.μ * sum(z[r] for r  = 1:initial_parameters.num_const )
            )

            +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[:], s ) .* x )
            #+  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[:], s ) * x[i] for i = 1:initial_parameters.num_first_stage_var)

        )

        #quadratic constraints
        @constraint( vector_of_subproblems[s], [ r = 1 : initial_parameters.num_const ],
            sum( y[i] * generated_parameters.constraint_Qs[s][r][i, j] * y[j] for i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var)
            + sum( x[i] * generated_parameters.constraint_fs[s][r][1, i] for i = 1:initial_parameters.num_first_stage_var)
            + sum( y[j] * generated_parameters.constraint_fs[s][r][2, j] for j = 1:initial_parameters.num_second_stage_var )
            + generated_parameters.constraint_fs[s][r][3, 1] + z[r] <= 0
        )

        # box constraints for first stage variables
        @constraint( vector_of_subproblems[s],
            generated_parameters.x_boundaries[:, 1] .<= x .<= generated_parameters.x_boundaries[:, 2])

        # box constraints for second stage variables
        @constraint( vector_of_subproblems[s],
            generated_parameters.y_boundaries[s][:, 1] .<= y .<= generated_parameters.y_boundaries[s][:, 2])

    end

    return vector_of_subproblems

end


function RNMDT_based_lagrangian_relaxation_problem_generation(initial_parameters::MIP_initial_parameters, generated_parameters::MIP_generated_parameters, precision_p::Array{Int64})

    # generate the starting values for the lagrangian multipliers
    vector_of_lambda_lagrangian = Array{Array{Float64}}(undef, initial_parameters.num_scen - 1)
    [ vector_of_lambda_lagrangian[i] = zeros(1,initial_parameters.num_first_stage_var)
            for i = 1 : initial_parameters.num_scen - 1 ]

    # creating the array of subproblems
    vector_of_subproblems = Array{Any}(undef, initial_parameters.num_scen)

    # formulating the subproblems based on the scenarios
    for s = 1 : initial_parameters.num_scen

        vector_of_subproblems[s] = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol, "Method" => initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag,  "Threads" => initial_parameters.gurobi_parameters.Threads, "NumericFocus" => initial_parameters.gurobi_parameters.NumericFocus))

        # first stage decision variables
        @variable( vector_of_subproblems[s], x[1 : initial_parameters.num_first_stage_var], Int )
        JuMP.unset_integer.(vector_of_subproblems[s][:x][generated_parameters.x_cont_indexes])

        # second stage decision variables
        @variable(vector_of_subproblems[s], y[1 : initial_parameters.num_second_stage_var], Int )
        JuMP.unset_integer.(vector_of_subproblems[s][:y][generated_parameters.y_cont_indexes[:,s]])

        # slack variables
        @variable(vector_of_subproblems[s], z[ 1 : initial_parameters.num_const ] >=0 )

        # RNMDT variables
        @variables vector_of_subproblems[s] begin
            y_heat_RNMDT[ 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_second_stage_var, minimum(precision_p) : -1 ] # x_heat
            delta_y_RNMDT[ 1 : initial_parameters.num_second_stage_var ] # delta_y_heat
            z_RNMDT[ 1 : initial_parameters.num_second_stage_var, minimum(precision_p) : -1 ], Bin
            w_RNMDT[ 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_second_stage_var ]
            delta_w_RNMDT[1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_second_stage_var ]

        end

        # quadratic objective
        @objective( vector_of_subproblems[s], Min,
            - initial_parameters.scen_prob[s] *
            ( sum(generated_parameters.objective_Qs[s][i, j] * w_RNMDT[i, j]
                for i = 1 : initial_parameters.num_second_stage_var,
                    j = 1 : initial_parameters.num_second_stage_var)
            + sum( x[i] * generated_parameters.objective_c[i]  for i = 1:initial_parameters.num_first_stage_var)
            + sum( y[j] * generated_parameters.objective_fs[s][j]  for j = 1:initial_parameters.num_second_stage_var)
            - initial_parameters.μ * sum(z[r] for r  = 1:initial_parameters.num_const )
            )

            +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[:], s ) .* x )

        )

        #quadratic constraint
        @constraint( vector_of_subproblems[s], [ r = 1 : initial_parameters.num_const ],
            sum(generated_parameters.constraint_Qs[s][r][i, j] * w_RNMDT[i, j]
                for i = 1 : initial_parameters.num_second_stage_var,
                j = 1 : initial_parameters.num_second_stage_var)
            + sum( x[i] * generated_parameters.constraint_fs[s][r][1, i] for i = 1:initial_parameters.num_first_stage_var)
            + sum( y[j] * generated_parameters.constraint_fs[s][r][2, j] for j = 1:initial_parameters.num_second_stage_var )
            + generated_parameters.constraint_fs[s][r][3, 1] + z[r] <= 0
        ) # 27

        # box constraints for first stage variables
        @constraint( vector_of_subproblems[s],
            generated_parameters.x_boundaries[:, 1] .<= x .<= generated_parameters.x_boundaries[:, 2])

        # box constraints for second stage variables
        @constraint( vector_of_subproblems[s],
            generated_parameters.y_boundaries[s][:, 1] .<= y .<= generated_parameters.y_boundaries[s][:, 2])

## auxiliary RNMDT constraints
        @constraint( vector_of_subproblems[s], [ j  = 1 : initial_parameters.num_second_stage_var],
            y[j] == (generated_parameters.y_boundaries[s][j, 2] - generated_parameters.y_boundaries[s][j, 1])  *
            ( sum( 2.0^l * z_RNMDT[j, l] for l = precision_p[j, s] : -1) + delta_y_RNMDT[j] ) ) # 28

        @constraint( vector_of_subproblems[s], [ i  = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var ],
            w_RNMDT[i, j] == (generated_parameters.y_boundaries[s][j, 2] - generated_parameters.y_boundaries[s][j, 1]) *
            ( sum( 2.0^l * y_heat_RNMDT[i, j, l] for l = precision_p[j, s] : -1) + delta_w_RNMDT[i, j] ) ) # 29

        @constraint( vector_of_subproblems[s], [ j  = 1 : initial_parameters.num_second_stage_var ],
            0 <= delta_y_RNMDT[j] <= 2.0^precision_p[j, s] ) # 30

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var],
            2.0^precision_p[i, s] * ( y[i] - generated_parameters.y_boundaries[s][i, 2] ) + generated_parameters.y_boundaries[s][i, 2] * delta_y_RNMDT[j] <= delta_w_RNMDT[i, j]) # 31

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var],
            delta_w_RNMDT[i, j] <= 2.0^precision_p[i, s] * ( y[i] - generated_parameters.y_boundaries[s][i, 1] ) + generated_parameters.y_boundaries[s][i, 1] * delta_y_RNMDT[j]) # 31

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var],
            generated_parameters.y_boundaries[s][i, 1] * delta_y_RNMDT[j] <= delta_w_RNMDT[i, j]) #32

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var],
            delta_w_RNMDT[i, j] <= generated_parameters.y_boundaries[s][i, 2] * delta_y_RNMDT[j]) # 32

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, l = precision_p[j, s] : -1],
            generated_parameters.y_boundaries[s][i, 1] * z_RNMDT[j, l]  <= y_heat_RNMDT[i, j, l]) # 33

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, l = precision_p[j, s] : -1],
            y_heat_RNMDT[i, j, l] <= generated_parameters.y_boundaries[s][i, 2] * z_RNMDT[j, l]) # 33

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, l = precision_p[j, s] : -1],
            generated_parameters.y_boundaries[s][i, 1] * (1 - z_RNMDT[j, l]) <= y[i] - y_heat_RNMDT[i, j, l] ) # 34

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, l = precision_p[j, s] : -1],
            y[i] - y_heat_RNMDT[i, j, l] <= generated_parameters.y_boundaries[s][i, 2] * (1 - z_RNMDT[j, l]) ) # 34
##
    end

    return vector_of_subproblems

end

function RNMDT_based_augmented_lagrangian_relaxation_problem_generation(initial_parameters::MIP_initial_parameters, generated_parameters::MIP_generated_parameters, precision_p::Array{Int64})

    # generate the starting values for the lagrangian multipliers
    vector_of_lambda_lagrangian = Array{Array{Float64}}(undef, initial_parameters.num_scen)
    [ vector_of_lambda_lagrangian[i] = zeros(1,initial_parameters.num_first_stage_var)
            for i = 1 : initial_parameters.num_scen ]

    # creating the array of subproblems
    vector_of_subproblems = Array{Any}(undef, initial_parameters.num_scen)

    # formulating the subproblems based on the scenarios
    for s = 1 : initial_parameters.num_scen

        vector_of_subproblems[s] = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol, "Method" => initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag,  "Threads" => initial_parameters.gurobi_parameters.Threads, "NumericFocus" => initial_parameters.gurobi_parameters.NumericFocus))

        # first stage decision variables
        @variable(vector_of_subproblems[s], x[1 : initial_parameters.num_first_stage_var], Int )
        JuMP.unset_integer.(vector_of_subproblems[s][:x][generated_parameters.x_cont_indexes])

        # non-anticipativity conditions variable (augmented lagrangian)
        @variable(vector_of_subproblems[s], al_z[1 : initial_parameters.num_first_stage_var], Int )
        JuMP.unset_integer.(vector_of_subproblems[s][:al_z][generated_parameters.x_cont_indexes])

        # second stage decision variables
        @variable(vector_of_subproblems[s], y[1 : initial_parameters.num_second_stage_var], Int )
        JuMP.unset_integer.(vector_of_subproblems[s][:y][generated_parameters.y_cont_indexes[:,s]])

        # slack variables
        @variable(vector_of_subproblems[s], z[ 1 : initial_parameters.num_const ] >=0 )

        # RNMDT variables
        @variables vector_of_subproblems[s] begin
            y_heat_RNMDT[ 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_second_stage_var, minimum(precision_p) : -1 ] # x_heat
            delta_y_RNMDT[ 1 : initial_parameters.num_second_stage_var ] # delta_y_heat
            z_RNMDT[ 1 : initial_parameters.num_second_stage_var, minimum(precision_p) : -1 ], Bin
            w_RNMDT[ 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_second_stage_var ]
            delta_w_RNMDT[1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_second_stage_var ]

        end

        # quadratic objective
        @objective( vector_of_subproblems[s], Min,
            - initial_parameters.scen_prob[s] *
            ( sum(generated_parameters.objective_Qs[s][i, j] * w_RNMDT[i, j]
                for i = 1 : initial_parameters.num_second_stage_var,
                    j = 1 : initial_parameters.num_second_stage_var)
            + sum( x[i] * generated_parameters.objective_c[i]  for i = 1:initial_parameters.num_first_stage_var)
            + sum( y[j] * generated_parameters.objective_fs[s][j]  for j = 1:initial_parameters.num_second_stage_var)
            - initial_parameters.μ * sum(z[r] for r  = 1:initial_parameters.num_const )
            )
            + sum( vector_of_lambda_lagrangian[s] .* (x .- al_z) )

            + sum( initial_parameters.al_penalty_parameter/2 .* (x .- al_z) .* (x .- al_z) )
            #+  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[:], s ) .* x )

        )

        #quadratic constraint
        @constraint( vector_of_subproblems[s], [ r = 1 : initial_parameters.num_const ],
            sum(generated_parameters.constraint_Qs[s][r][i, j] * w_RNMDT[i, j]
                for i = 1 : initial_parameters.num_second_stage_var,
                j = 1 : initial_parameters.num_second_stage_var)
            + sum( x[i] * generated_parameters.constraint_fs[s][r][1, i] for i = 1:initial_parameters.num_first_stage_var)
            + sum( y[j] * generated_parameters.constraint_fs[s][r][2, j] for j = 1:initial_parameters.num_second_stage_var )
            + generated_parameters.constraint_fs[s][r][3, 1] - z[r] <= 0
        ) # 27

        # box constraints for first stage variables
        @constraint( vector_of_subproblems[s],
            generated_parameters.x_boundaries[:, 1] .<= x .<= generated_parameters.x_boundaries[:, 2])

        # box constraints for second stage variables
        @constraint( vector_of_subproblems[s],
            generated_parameters.y_boundaries[s][:, 1] .<= y .<= generated_parameters.y_boundaries[s][:, 2])

## auxiliary RNMDT constraints
        @constraint( vector_of_subproblems[s], [ j  = 1 : initial_parameters.num_second_stage_var],
            y[j] == (generated_parameters.y_boundaries[s][j, 2] - generated_parameters.y_boundaries[s][j, 1])  *
            ( sum( 2.0^l * z_RNMDT[j, l] for l = precision_p[j, s] : -1) + delta_y_RNMDT[j] ) ) # 28

        @constraint( vector_of_subproblems[s], [ i  = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var ],
            w_RNMDT[i, j] == (generated_parameters.y_boundaries[s][j, 2] - generated_parameters.y_boundaries[s][j, 1]) *
            ( sum( 2.0^l * y_heat_RNMDT[i, j, l] for l = precision_p[j, s] : -1) + delta_w_RNMDT[i, j] ) ) # 29

        @constraint( vector_of_subproblems[s], [ j  = 1 : initial_parameters.num_second_stage_var ],
            0 <= delta_y_RNMDT[j] <= 2.0^precision_p[j, s] ) # 30

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var],
            2.0^precision_p[i, s] * ( y[i] - generated_parameters.y_boundaries[s][i, 2] ) + generated_parameters.y_boundaries[s][i, 2] * delta_y_RNMDT[j] <= delta_w_RNMDT[i, j]) # 31

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var],
            delta_w_RNMDT[i, j] <= 2.0^precision_p[i, s] * ( y[i] - generated_parameters.y_boundaries[s][i, 1] ) + generated_parameters.y_boundaries[s][i, 1] * delta_y_RNMDT[j]) # 31

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var],
            generated_parameters.y_boundaries[s][i, 1] * delta_y_RNMDT[j] <= delta_w_RNMDT[i, j]) #32

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var],
            delta_w_RNMDT[i, j] <= generated_parameters.y_boundaries[s][i, 2] * delta_y_RNMDT[j]) # 32

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, l = precision_p[j, s] : -1],
            generated_parameters.y_boundaries[s][i, 1] * z_RNMDT[j, l]  <= y_heat_RNMDT[i, j, l]) # 33

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, l = precision_p[j, s] : -1],
            y_heat_RNMDT[i, j, l] <= generated_parameters.y_boundaries[s][i, 2] * z_RNMDT[j, l]) # 33

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, l = precision_p[j, s] : -1],
            generated_parameters.y_boundaries[s][i, 1] * (1 - z_RNMDT[j, l]) <= y[i] - y_heat_RNMDT[i, j, l] ) # 34

        @constraint( vector_of_subproblems[s], [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, l = precision_p[j, s] : -1],
            y[i] - y_heat_RNMDT[i, j, l] <= generated_parameters.y_boundaries[s][i, 2] * (1 - z_RNMDT[j, l]) ) # 34
##
    end

    return vector_of_subproblems

end

function RNMDT_based_problem_generation(initial_parameters::MIP_initial_parameters, generated_parameters::MIP_generated_parameters)

    # all the commments for the auxiliary constraints in this section refer to the paper
    # Enhancing the normalized multiparametric disaggregation
    # technique for mixed-integer quadratic programming

    # simplifying the notations
    precision_p = initial_parameters.RNMDT_precision_factor

    # defining RNMDT problem
    RNMDT_problem = Model( optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol, "Method" => initial_parameters.gurobi_parameters.Method, "Threads" => initial_parameters.gurobi_parameters.Threads, "TimeLimit" => initial_parameters.solver_time_limit, "Presolve" => 0) )

    # first stage decision variables
    @variable(RNMDT_problem, x[ 1 : initial_parameters.num_first_stage_var, 1 : initial_parameters.num_scen ], Int )
    JuMP.unset_integer.(RNMDT_problem[:x][generated_parameters.x_cont_indexes, :])

    # second stage decision variables
    @variable(RNMDT_problem, y[ 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_scen ], Int)
    [JuMP.unset_integer.(RNMDT_problem[:y][generated_parameters.y_cont_indexes[:,j], j]) for j = 1:initial_parameters.num_scen]

    # slack variables
    @variable(RNMDT_problem, z[ 1 : initial_parameters.num_scen, 1 : initial_parameters.num_const ] >=0 )

    #auxiliary RNMDT variables
    @variables RNMDT_problem begin
        y_heat_RNMDT[ 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_second_stage_var, minimum(precision_p) : -1, 1 : initial_parameters.num_scen ] # x_heat
        delta_y_RNMDT[ 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_scen ] # delta_x_heat
        z_RNMDT[ 1 : initial_parameters.num_second_stage_var, minimum(precision_p) : -1, 1 : initial_parameters.num_scen ], Bin
        w_RNMDT[ 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_scen ]
        delta_w_RNMDT[1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_scen ]
    end

    #quadratic objective
    @objective( RNMDT_problem, Min,
        - sum( initial_parameters.scen_prob[s] *
            (
            sum(generated_parameters.objective_Qs[s][i, j] * w_RNMDT[i, j, s]
                    for i = 1 : initial_parameters.num_second_stage_var,
                        j = 1 : initial_parameters.num_second_stage_var)
                + sum( x[i, s] * generated_parameters.objective_c[i] for i = 1:initial_parameters.num_first_stage_var)
                + sum( y[j, s] * generated_parameters.objective_fs[s][j] for j = 1:initial_parameters.num_second_stage_var)
                - initial_parameters.μ * sum(z[s,r] for r  = 1:initial_parameters.num_const )
            )
        for s = 1 : initial_parameters.num_scen)
    ) # 26

    #quadratic constraints
    @constraint( RNMDT_problem, [ s = 1 : initial_parameters.num_scen, r = 1 : initial_parameters.num_const ],
        sum(generated_parameters.constraint_Qs[s][r][i, j] * w_RNMDT[i, j, s]
            for i = 1 : initial_parameters.num_second_stage_var,
                j = 1 : initial_parameters.num_second_stage_var)
        + sum( x[i, s] * generated_parameters.constraint_fs[s][r][1, i] for i = 1:initial_parameters.num_first_stage_var)
        + sum( y[j, s] * generated_parameters.constraint_fs[s][r][2, j] for j = 1:initial_parameters.num_second_stage_var )
        + generated_parameters.constraint_fs[s][r][3, 1] + z[s,r] <= 0
    ) # 27

## auxiliary RNMDT constraints

    @constraint( RNMDT_problem, [ j  = 1 : initial_parameters.num_second_stage_var, s  = 1 : initial_parameters.num_scen],
        y[j, s] == (generated_parameters.y_boundaries[s][j, 2] - generated_parameters.y_boundaries[s][j, 1])  *
        ( sum( 2.0^l * z_RNMDT[j, l, s] for l = precision_p[j, s] : -1) + delta_y_RNMDT[j, s] ) ) # 28

    @constraint( RNMDT_problem, [ i  = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen ],
        w_RNMDT[i, j, s] == (generated_parameters.y_boundaries[s][j, 2] - generated_parameters.y_boundaries[s][j, 1]) *
        ( sum( 2.0^l * y_heat_RNMDT[i, j, l, s] for l = precision_p[j, s] : -1) + delta_w_RNMDT[i, j, s] ) ) # 29

    @constraint( RNMDT_problem, [ j  = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen ],
        0 <= delta_y_RNMDT[j, s] <= 2.0^precision_p[j, s] ) # 30

    @constraint( RNMDT_problem, [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen],
        2.0^precision_p[i, s] * ( y[i, s] - generated_parameters.y_boundaries[s][i, 2] ) + generated_parameters.y_boundaries[s][i, 2] * delta_y_RNMDT[j, s] <= delta_w_RNMDT[i, j, s]) # 31

    @constraint( RNMDT_problem, [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen],
        delta_w_RNMDT[i, j, s] <= 2.0^precision_p[i, s] * ( y[i, s] - generated_parameters.y_boundaries[s][i, 1] ) + generated_parameters.y_boundaries[s][i, 1] * delta_y_RNMDT[j, s]) # 31

    @constraint( RNMDT_problem, [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen ],
        generated_parameters.y_boundaries[s][i, 1] * delta_y_RNMDT[j, s] <= delta_w_RNMDT[i, j, s]) #32

    @constraint( RNMDT_problem, [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen],
        delta_w_RNMDT[i, j, s] <= generated_parameters.y_boundaries[s][i, 2] * delta_y_RNMDT[j, s]) # 32

    @constraint( RNMDT_problem, [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen, l = precision_p[j, s] : -1 ],
        generated_parameters.y_boundaries[s][i, 1] * z_RNMDT[j, l, s]  <= y_heat_RNMDT[i, j, l, s]) # 33

    @constraint( RNMDT_problem, [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen, l = precision_p[j, s] : -1 ],
        y_heat_RNMDT[i, j, l, s] <= generated_parameters.y_boundaries[s][i, 2] * z_RNMDT[j, l, s]) # 33

    @constraint( RNMDT_problem, [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen, l = precision_p[j, s] : -1 ],
        generated_parameters.y_boundaries[s][i, 1] * (1 - z_RNMDT[j, l, s]) <= y[i,s] - y_heat_RNMDT[i, j, l, s] ) # 34

    @constraint( RNMDT_problem, [ i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen, l = precision_p[j, s] : -1 ],
        y[i,s] - y_heat_RNMDT[i, j, l, s] <= generated_parameters.y_boundaries[s][i, 2] * (1 - z_RNMDT[j, l, s]) ) # 34

##
    # box constraints for first stage variables
    @constraint(RNMDT_problem, [ s = 1 : initial_parameters.num_scen ],
        generated_parameters.x_boundaries[:, 1] .<= x[:, s] .<= generated_parameters.x_boundaries[:, 2])

    # box constraints for second stage variables
    @constraint(RNMDT_problem, [s = 1 : initial_parameters.num_scen],
        generated_parameters.y_boundaries[s][:, 1] .<= y[:, s] .<= generated_parameters.y_boundaries[s][:, 2])

    # non-anticipativity conditions
    @constraint(RNMDT_problem, [s in 2 : initial_parameters.num_scen, i = 1 : initial_parameters.num_first_stage_var],
        x[i, s] - x[i, 1] == 0 )

    return RNMDT_problem
end


function primal_problem_based_McCormic_relaxation(initial_parameters::MIP_initial_parameters, generated_parameters::MIP_generated_parameters)

        # defining RNMDT problem
        McCormick_relaxation = Model( optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol, "Method" => initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag, "Threads" => initial_parameters.gurobi_parameters.Threads, "TimeLimit" => initial_parameters.solver_time_limit) )

        # first stage decision variables
        @variable(McCormick_relaxation, x[ 1 : initial_parameters.num_first_stage_var, 1 : initial_parameters.num_scen ], Int )
        JuMP.unset_integer.(McCormick_relaxation[:x][generated_parameters.x_cont_indexes, :])

        # second stage decision variables
        @variable(McCormick_relaxation, y[ 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_scen ], Int)
        [JuMP.unset_integer.(McCormick_relaxation[:y][generated_parameters.y_cont_indexes[:,j], j]) for j = 1:initial_parameters.num_scen]

        # slack variables
        @variable(McCormick_relaxation, z[ 1 : initial_parameters.num_scen, 1 : initial_parameters.num_const ] >=0 )

        #auxiliary McCormicj relaxation variables
        @variables McCormick_relaxation begin
            w_McCormick[ 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_scen ]
        end

        #quadratic objective
        @objective( McCormick_relaxation, Min,
            - sum( initial_parameters.scen_prob[s] *
                (
                sum(generated_parameters.objective_Qs[s][i, j] * w_McCormick[i, j, s]
                        for i = 1 : initial_parameters.num_second_stage_var,
                            j = 1 : initial_parameters.num_second_stage_var)
                    + sum( x[i, s] * generated_parameters.objective_c[i] for i = 1:initial_parameters.num_first_stage_var)
                    + sum( y[j, s] * generated_parameters.objective_fs[s][j] for j = 1:initial_parameters.num_second_stage_var)
                    - initial_parameters.μ * sum(z[s,r] for r  = 1:initial_parameters.num_const )
                )
            for s = 1 : initial_parameters.num_scen)
        ) # 26

        #quadratic constraints
        @constraint( McCormick_relaxation, [ s = 1 : initial_parameters.num_scen, r = 1 : initial_parameters.num_const ],
            sum(generated_parameters.constraint_Qs[s][r][i, j] * w_McCormick[i, j, s]
                for i = 1 : initial_parameters.num_second_stage_var,
                    j = 1 : initial_parameters.num_second_stage_var)
            + sum( x[i, s] * generated_parameters.constraint_fs[s][r][1, i] for i = 1:initial_parameters.num_first_stage_var)
            + sum( y[j, s] * generated_parameters.constraint_fs[s][r][2, j] for j = 1:initial_parameters.num_second_stage_var )
            + generated_parameters.constraint_fs[s][r][3, 1] + z[s,r] <= 0
        ) # 27

    ## auxiliary McCormick relaxation constraints

        @constraint( McCormick_relaxation, [ i  = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen ],
            w_McCormick[i, j, s] >= y[i, s]*generated_parameters.y_boundaries[s][j, 1] + y[j, s]*generated_parameters.y_boundaries[s][i, 1] - generated_parameters.y_boundaries[s][i, 1]*generated_parameters.y_boundaries[s][j, 1])

        @constraint( McCormick_relaxation, [ i  = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen ],
            w_McCormick[i, j, s] >= y[i, s]*generated_parameters.y_boundaries[s][j, 2] + y[j, s]*generated_parameters.y_boundaries[s][i, 2] - generated_parameters.y_boundaries[s][i, 2]*generated_parameters.y_boundaries[s][j, 2])

        @constraint( McCormick_relaxation, [ i  = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen ],
            w_McCormick[i, j, s] <= y[i, s]*generated_parameters.y_boundaries[s][j, 1] + y[j, s]*generated_parameters.y_boundaries[s][i, 2] - generated_parameters.y_boundaries[s][i, 2]*generated_parameters.y_boundaries[s][j, 1])

        @constraint( McCormick_relaxation, [ i  = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen ],
            w_McCormick[i, j, s] <= y[i, s]*generated_parameters.y_boundaries[s][j, 2] + y[j, s]*generated_parameters.y_boundaries[s][i, 1] - generated_parameters.y_boundaries[s][i, 1]*generated_parameters.y_boundaries[s][j, 2])

    ##
        # box constraints for first stage variables
        @constraint(McCormick_relaxation, [ s = 1 : initial_parameters.num_scen ],
            generated_parameters.x_boundaries[:, 1] .<= x[:, s] .<= generated_parameters.x_boundaries[:, 2])

        # box constraints for second stage variables
        @constraint(McCormick_relaxation, [s = 1 : initial_parameters.num_scen],
            generated_parameters.y_boundaries[s][:, 1] .<= y[:, s] .<= generated_parameters.y_boundaries[s][:, 2])

        # non-anticipativity conditions
        @constraint(McCormick_relaxation, [s in 2 : initial_parameters.num_scen, i = 1 : initial_parameters.num_first_stage_var],
            x[i, s] - x[i, 1] == 0 )

    return McCormick_relaxation
end
