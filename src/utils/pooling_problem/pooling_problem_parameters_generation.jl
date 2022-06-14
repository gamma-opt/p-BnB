## This file contains the fucntion required to generate parameters for the pooling problem generation 
## in a form suitable for the input in the solver

function polling_problem_parameters_generation( model::JuMP.Model, number_of_scenarios::Int, flow_cost, pool_cost, precision_p, bub)
    Variables, Objective_c, Affine_constraints, Quadratic_constraints, Box_constraints, bub = polling_problem_parameters_extraction(model, bub)

    second_stage_variables_number = length(Variables)

    # gathering the variables inidices that correpspond to all arcs
    arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) ,Variables)]
    
    #saving indices
    arc_var_ind = Array{Int64}(undef,length(arc_var))
    [arc_var_ind[i] = arc_var[i].index for i = 1:length(arc_var)]

    number_of_pools = length(findall(x -> occursin("q[p", x.name), Variables))

    # then the number of frist-stage variables would be the number of pools + the number of arcs
    fs_num = number_of_pools + length(arc_var)
    #@show fs_num

    # generating the MIP_initital_parameters structure  
    initial_parameters = initialisation(number_of_scenarios, fs_num, length(Variables),length(Quadratic_constraints), precision_p)
            
    #---------------generating integer and continuous indexes-----------------------
    Random.seed!(initial_parameters.random_seed) 

    shuffled_x_ind = randperm!(Random.seed!(initial_parameters.random_seed), collect(1:initial_parameters.num_first_stage_var))
    x_int_indexes = shuffled_x_ind[1:initial_parameters.num_int_var_first_stage]
    x_cont_indexes = shuffled_x_ind[initial_parameters.num_int_var_first_stage+1:initial_parameters.num_first_stage_var]
    
    shuffled_y_ind = Array{Int,2}(undef,initial_parameters.num_second_stage_var, initial_parameters.num_scen )
    #@show randperm!(Random.seed!(initial_parameters.random_seed + 1 - 1), collect(1:initial_parameters.num_second_stage_var))
    [shuffled_y_ind[:, j] = randperm!(Random.seed!(initial_parameters.random_seed + j - 1), collect(1:initial_parameters.num_second_stage_var)) for j = 1:initial_parameters.num_scen]
    y_int_indexes = shuffled_y_ind[1:initial_parameters.num_int_var_second_stage, :]
     y_cont_indexes = shuffled_y_ind[initial_parameters.num_int_var_second_stage+1:initial_parameters.num_second_stage_var, :]
    

    #---------------generating box constraints for the scenarios--------------

    # box constraints for the second-stage Variables
    # At first we generate them for a single scenario
    
    ss_box_constraint_left = Array{Float64}(undef, length(Variables))
    ss_box_constraint_right = Array{Float64}(undef, length(Variables))

    # Extracting existing left hand side values 
    for i = 1:length(Box_constraints)
        ss_box_constraint_left[Box_constraints[i].variable_index] = Box_constraints[i].RHS
    end

    # splitting affine constraints in equality and inequality constraints
    Affine_equality_constraints = Affine_constraints[findfirst(x->occursin("EqualTo", x.sign), Affine_constraints)]
    Affine_inequality_constraints = Affine_constraints[findfirst(x->occursin("LessThan", x.sign), Affine_constraints)]

    # Extracting existing right hand side values from affine ineqaulity constraints 
    for i = 1:size(Affine_inequality_constraints.linear_coefficients_matrix,1)
        # if we have the constraint with only one variable
        if length(findall(x->x!=0,Affine_inequality_constraints.linear_coefficients_matrix[i,:])) == 1 
            ss_box_constraint_right[findfirst(x -> x!= 0, Affine_inequality_constraints.linear_coefficients_matrix[i,:])] = Affine_inequality_constraints.RHS[i]
        else
            # collect all the non-zero indices 
            indices =  findall(x -> x!= 0, Affine_inequality_constraints.linear_coefficients_matrix[i,:])
            ss_box_constraint_right[indices] .= Affine_inequality_constraints.RHS[i]
        end
    end


    # Extracting right hand side values from affine equality constraints 
    for i = 1:size(Affine_equality_constraints.linear_coefficients_matrix,1)
        # if we have the constraint with only one variable
        if length(findall(x->x!=0,Affine_equality_constraints.linear_coefficients_matrix[i,:])) == 1 
            ss_box_constraint_right[findfirst(x -> x!= 0, Affine_equality_constraints.linear_coefficients_matrix[i,:])] = Affine_equality_constraints.RHS[i]
            ss_box_constraint_left[findfirst(x -> x!= 0, Affine_equality_constraints.linear_coefficients_matrix[i,:])] = Affine_equality_constraints.RHS[i]
        else 
            # collect all negative indicies 
            negative_indices =  findall(x -> x < 0, Affine_equality_constraints.linear_coefficients_matrix[i,:])
            # collect all positive indicies 
            postive_indices =  findall(x -> x > 0, Affine_equality_constraints.linear_coefficients_matrix[i,:])
            ss_box_constraint_right[postive_indices] .= sum(ss_box_constraint_right[negative_indices]) + Affine_equality_constraints.RHS[i]
        end

    end
    
    # Set the boudns of the quality value of the pooling node 

    # extract all the variables correspoding to quality values at pooling nodes
    qp_variables = Variables[findall(x -> occursin("q[p", x.name), Variables)]
        
    # extract all the variables correspoding to quality values 
    q_not_p_variables = Variables[findall(x -> occursin("q[", x.name)&&!occursin("q[p", x.name) , Variables)]

    # create an array of corresponding indices of the variables correspoding to quality values 
    q_not_p_variables_indices = Array{Int64}(undef,length(q_not_p_variables))
    [q_not_p_variables_indices[i] = q_not_p_variables[i].index for i = 1:length(q_not_p_variables)]

    # set the boudns as the max of all quality bounds and the min of all quality bounds respectively 
    for variable in qp_variables
        ss_box_constraint_right[variable.index] = maximum(ss_box_constraint_right[ q_not_p_variables_indices])
        ss_box_constraint_left[variable.index] = minimum(ss_box_constraint_left[ q_not_p_variables_indices])
    end

    # generating bounds for the second stage variables for all the scenarios
    y_boundaries = Array{Any}(undef, initial_parameters.num_scen)
    [y_boundaries[j] = [ss_box_constraint_left ss_box_constraint_right] for j = 1:initial_parameters.num_scen]
    
    # generating bounds for the first stage variables 
    x_boundaries = [0*ones( initial_parameters.num_first_stage_var, 1 ) 1*ones( initial_parameters.num_first_stage_var, 1 )]

    #---------------generating quadratic constraints for the scenarios--------------

    # generating matrices Qsi for the left hand side of the contraint for each of the scenario
    constraint_Qs = Array{Any}(undef, 1, initial_parameters.num_scen)

    [ constraint_Qs[i] = [0.5 .* (Quadratic_constraints[j].quadratic_coefficients_matrix + Quadratic_constraints[j].quadratic_coefficients_matrix') 
        for j = 1 : initial_parameters.num_const] for i = 1:initial_parameters.num_scen ]

    # generating affine functions' coefficients for the left hand side of the constraint for each of the scenario
    constraint_fs = Array{Any}(undef, 1, initial_parameters.num_scen)

    [ constraint_fs[i] = [  [  [zeros(1, initial_parameters.num_first_stage_var) zeros(1, ( initial_parameters.num_second_stage_var > initial_parameters.num_first_stage_var ) ? (initial_parameters.num_second_stage_var - initial_parameters.num_first_stage_var) : 0 ) ];
                                [ Quadratic_constraints[j].affine_coefficients_vector' zeros(1, ( initial_parameters.num_first_stage_var > initial_parameters.num_second_stage_var ) ? (initial_parameters.num_first_stage_var - initial_parameters.num_second_stage_var) : 0 ) ];
                                Quadratic_constraints[j].constant.* [1 zeros(1, ( initial_parameters.num_second_stage_var > initial_parameters.num_first_stage_var ) ? (initial_parameters.num_second_stage_var - 1) : (initial_parameters.num_first_stage_var - 1))] 
                            ]
                         for j = 1:initial_parameters.num_const] for 
    i = 1:initial_parameters.num_scen ]
    # first row - x_coeficients (integer variables)
    # second row - y_coeficients (continuous variables)
    # third row  - constant

    #---------------generating obejctive fucntions for the scenarios----------------

    # generating matrices Qsi for the objecyive for each of the scenario
    objective_Qs = Array{Any}(undef, 1, initial_parameters.num_scen)

    [ objective_Qs[i] = -1. * zeros(initial_parameters.num_second_stage_var, initial_parameters.num_second_stage_var)
        for i = 1:initial_parameters.num_scen ]

    # generating linear functions' coefficients for the objective for each of the scenario
    objective_fs = Array{Any}(undef, 1, initial_parameters.num_scen)

    #@show Objective_c

    # for the second stage variables 
    [ objective_fs[i] = -1 .* Objective_c for i = 1:initial_parameters.num_scen  ]
        
    # for the first-stage variables (scenario-independent coefficients)
    objective_c = -1 .* [-flow_cost .* ones(1, length(arc_var)) -pool_cost .* ones(1, number_of_pools)]
    # first arc related costs then pool related costs

    generated_parameters = MIP_generated_parameters(constraint_Qs, constraint_fs, objective_Qs, objective_fs, objective_c, x_boundaries, x_int_indexes, x_cont_indexes, y_boundaries, y_int_indexes, y_cont_indexes)

    # if we need to print out the correlation map between pooling problem variables and resulting problem variables indicies
    if initial_parameters.pool_problem_is_used && initial_parameters.pool_prob_par.print_out_variables_correlation_map
        # create an array that will contain left part 
        pp_index = []
        # create an array that will contain right part 
        stoch_pp_index = []

        push!(pp_index, "second-stage")
        push!(stoch_pp_index, "variables")

        for i = 1:length(Variables)
            if Variables[i] in arc_var
                push!(pp_index, "flow: " * Variables[i].name)
            else 
                push!(pp_index, Variables[i].name)
            end
            push!(stoch_pp_index, "y["*string(Variables[i].index)*"]")
        end
        push!(pp_index, "first-stage")
        push!(stoch_pp_index, "binary variables")

        for i = 1:length(arc_var)
            push!(pp_index, "decision on: " * arc_var[i].name)
            push!(stoch_pp_index, "x["*string(i)*"]")
        end
        
        for i = 1:number_of_pools
            push!(pp_index, "decision on: " * "p"*string(i))
            push!(stoch_pp_index, "x["*string(length(arc_var)+i)*"]")
        end

        map_df = DataFrame(pooling_problem=pp_index, stoch_pooling_problem=stoch_pp_index)
        XLSX.writetable(initial_parameters.pool_prob_par.primal_pool_problem_link*"/map_var_name_index_"*string(Dates.today())*".xlsx",map_df, overwrite=true)
    end 

    return  initial_parameters, generated_parameters, Variables, bub
end 