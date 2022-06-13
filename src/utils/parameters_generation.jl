
#--------------auxiliary function for generating quadratic matrices for cosnstraints and obejctive with predefined densiity
function quadratic_matrix_generation(density, dimention, min_range, max_range, PSD, seed)

    # if the matrix is supposed to be PSD
    if PSD == "yes"

        # creating the matrix of predefined density,
        # taking into account that this density is defined for the case when a matrix is symmetrical
        # we generate triangular matrix in such a way that if it is represented as a symmetric one
        # (dividing upper triangular part by 2 and reflecting it towards the diagonal)
        # it will have that predefined density

        # separately generating diagonal elements not equal to zero (to be able to create PSD after)
        diagonal_elements  = max_range .* round.(rand(1, dimention), RoundUp, digits = 1)

        # calculating the number of non-zero elements in the upper diagonal part
        # depending on the predefined density
        # and taking into account that diagonal elements are non-zero
        number_of_other_non_zero_elements = max(0,  Int(round((density * dimention^2 - dimention) / 2)))

        # creating one-dimensional array containing an abovementioned number of non zero elements,
        # filling the rest with zeros and shuffling the array
        upper_triangular_elements  = shuffle!([(min_range .+ (max_range - min_range) .* round.(rand(1, number_of_other_non_zero_elements), digits = 1 )) zeros(1, Int(dimention*(dimention-1)/2) - number_of_other_non_zero_elements) ] )

        # creating final matrix with zero elements
        Q = zeros(dimention, dimention)

        # using above mentioned diagonal elements and upper triangular part
        # rewriting zeros in the upper triangular part accordingly
        for i = 1:dimention
            for j = i:dimention
                if i == j # if the element is on the diagonal using an array of diagonal elements
                    Q[i,i] = diagonal_elements[i]
                else
                    Q[i,j] = upper_triangular_elements[1]/2 # otherwise using one elemnts from upper triangular elements
                    Q[j,i] = upper_triangular_elements[1]/2 # otherwise using one elemnts from upper triangular elements
                    upper_triangular_elements = upper_triangular_elements[2:end] # deleting that element from the set
                end
            end
        end

        # if the matrix is supposed to be a PSD then check if the diagonal elements are bigger
        # than the sum of the elements on the same row and if it's not the case increase it by
        # a bigger enough number
        for i = 1:dimention
            #Q[i,i] = Q[i,i] > sum(Q[i,i+1:end]) ? Q[i,i] : sum(Q[i,i+1:end]) + 1
            Q[i,i] = Q[i,i] > (sum(abs.(Q[i, 1:i-1])) + sum(abs.(Q[i, i+1:end]))) ? Q[i,i] : (sum(abs.(Q[i, 1:i-1])) + sum(abs.(Q[i, i+1:end]))) +1000
        end
    # if the matrix does not have to be necessarily PSD
    else
        Q = 1 .* round.((min_range  .+ (max_range - min_range) .* Matrix(Symmetric(sprand(dimention,dimention,density)))), digits = 1)
    end

    return Q
end

#--------------function for generation the parameters---------------------------

function parameters_generation(initial_parameters::MIP_initial_parameters)

        Random.seed!(initial_parameters.random_seed)

        # defining max and min values for the quadratic matirces 
        # coefficients in constraints and objective
        #Max_value_for_quad_matrix_elements = 100
        #Min_value_for_quad_matrix_elements = 0


        # defining max and min values for first-stage variables 
        # linear coefficients in constraints and objective
        #Max_value_for_fs_linear_elements = 100
        #Min_value_for_fs_linear_elements = 0

        # defining max and min values for second-stage variables 
        # linear coefficients in constraints and objective
        #Max_value_for_ss_linear_elements = 100
        #Min_value_for_ss_linear_elements = 0
        

        # defining maximum value for the affine constant in the constraints
        #Min_value_for_affine_constant = 1000
        #Max_value_for_affine_constant = 100000

        #x_limits = [0 100] # max and min values for the continuous variables' box constraints
        #y_limits = [0 100] # max and min values for the integer variables' box constraints

        #---------------generating integer and continuous indexes-----------------------
        shuffled_x_ind = randperm!(Random.seed!(initial_parameters.random_seed), collect(1:initial_parameters.num_first_stage_var))
        x_int_indexes = shuffled_x_ind[1:initial_parameters.num_int_var_first_stage]
        x_cont_indexes = shuffled_x_ind[initial_parameters.num_int_var_first_stage+1:initial_parameters.num_first_stage_var]

        shuffled_y_ind = Array{Int,2}(undef,initial_parameters.num_second_stage_var, initial_parameters.num_scen )
        #@show randperm!(Random.seed!(initial_parameters.random_seed + 1 - 1), collect(1:initial_parameters.num_second_stage_var))
        [shuffled_y_ind[:, j] = randperm!(Random.seed!(initial_parameters.random_seed + j - 1), collect(1:initial_parameters.num_second_stage_var)) for j = 1:initial_parameters.num_scen]
        y_int_indexes = shuffled_y_ind[1:initial_parameters.num_int_var_second_stage, :]
        y_cont_indexes = shuffled_y_ind[initial_parameters.num_int_var_second_stage+1:initial_parameters.num_second_stage_var, :]

        #---------------generating constraints for the variables------------------------

        # in case x is Int (not Binary) we have a differnet interval defined for each x_i
        if !initial_parameters.bin_con_fs
            x_boundaries = [initial_parameters.box_con_fs_min*ones( initial_parameters.num_first_stage_var, 1 ) rand(
                Int( round( (initial_parameters.box_con_fs_max -initial_parameters.box_con_fs_min) / 2 ) ) :initial_parameters.box_con_fs_max,
                initial_parameters.num_first_stage_var, 1 ) ]
        else # In case x is Binary we have to consider [0,1] interval for all the variables x_i
            x_boundaries = [initial_parameters.box_con_fs_min*ones( initial_parameters.num_first_stage_var, 1 ) initial_parameters.box_con_fs_max*ones( initial_parameters.num_first_stage_var, 1 )]
        end 
        y_boundaries = Array{Any}(undef, initial_parameters.num_scen)
        [y_boundaries[j] = [initial_parameters.box_con_ss_min*ones( initial_parameters.num_second_stage_var, 1 ) rand(
             (initial_parameters.box_con_ss_max - initial_parameters.box_con_ss_min) / 2   : initial_parameters.box_con_ss_max,
            initial_parameters.num_second_stage_var, 1  ) ] for j = 1:initial_parameters.num_scen]

        #---------------generating quadratic constraints for the scenarios--------------

        # generating matrices Qsi for the left hand side of the contraint for each of the scenario
        constraint_Qs = Array{Any}(undef, 1, initial_parameters.num_scen)

        [ constraint_Qs[i] = [ quadratic_matrix_generation(initial_parameters.quad_mat_dens, initial_parameters.num_second_stage_var, initial_parameters.con_quad_min, initial_parameters.con_quad_max, "no", initial_parameters.random_seed)
            for j = 1 : initial_parameters.num_const] for i = 1:initial_parameters.num_scen ]

        # generating affine functions' coefficients for the left hand side of the constraint for each of the scenario
        constraint_fs = Array{Any}(undef, 1, initial_parameters.num_scen)

        [ constraint_fs[i] = [ 1 .* round.([(initial_parameters.con_lin_fs_min .+ (initial_parameters.con_lin_fs_max .- initial_parameters.con_lin_fs_min ) .* [rand(1, initial_parameters.num_first_stage_var) zeros(1, ( initial_parameters.num_second_stage_var > initial_parameters.num_first_stage_var ) ? (initial_parameters.num_second_stage_var - initial_parameters.num_first_stage_var) : 0 ) ]);
                                (initial_parameters.con_lin_ss_min .+ (initial_parameters.con_lin_ss_max - initial_parameters.con_lin_ss_min) .* [ rand(1, initial_parameters.num_second_stage_var) zeros(1, ( initial_parameters.num_first_stage_var > initial_parameters.num_second_stage_var ) ? (initial_parameters.num_first_stage_var - initial_parameters.num_second_stage_var) : 0 ) ]);
                                initial_parameters.con_aff_min .+ (initial_parameters.con_aff_max - initial_parameters.con_aff_min) .* [rand(1,1) zeros(1, ( initial_parameters.num_second_stage_var > initial_parameters.num_first_stage_var ) ? (initial_parameters.num_second_stage_var - 1) : (initial_parameters.num_first_stage_var - 1))] ], digits = 1)
            for j = 1:initial_parameters.num_const] for i = 1:initial_parameters.num_scen ]
        # first row - x_coeficients (integer variables)
        # second row - y_coeficients (continuous variables)
        # third row  - constant

        #---------------generating obejctive fucntions for the scenarios----------------

        # generating matrices Qsi for the objecyive for each of the scenario
        objective_Qs = Array{Any}(undef, 1, initial_parameters.num_scen)

        [ objective_Qs[i] = quadratic_matrix_generation(initial_parameters.quad_mat_dens, initial_parameters.num_second_stage_var, initial_parameters.obj_quad_min, initial_parameters.obj_quad_max, "no", initial_parameters.random_seed)
        for i = 1:initial_parameters.num_scen ]

        # generating linear functions' coefficients for the objective for each of the scenario
        objective_fs = Array{Any}(undef, 1, initial_parameters.num_scen)

        # for the second stage variables 
        [ objective_fs[i] = 1 .* round.( initial_parameters.obj_lin_ss_min .+ (initial_parameters.obj_lin_ss_max - initial_parameters.obj_lin_ss_min) .*  rand(1, initial_parameters.num_second_stage_var),  digits = 1) for i = 1:initial_parameters.num_scen  ]
        
        # for the first-stage variables (scenario-independent coefficients)
        objective_c = 1 .* round.( initial_parameters.obj_lin_fs_min .+ (initial_parameters.obj_lin_fs_min - initial_parameters.obj_lin_fs_min) .*  rand(1, initial_parameters.num_first_stage_var),  digits = 1)

        generated_parameters = MIP_generated_parameters(constraint_Qs, constraint_fs, objective_Qs, objective_fs, objective_c, x_boundaries, x_int_indexes, x_cont_indexes, y_boundaries, y_int_indexes, y_cont_indexes)

        return generated_parameters
end
