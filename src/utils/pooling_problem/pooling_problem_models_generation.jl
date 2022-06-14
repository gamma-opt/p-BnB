## This file contains all the necessary function required to deifine the modified stochastic 
## two-stage pooling problems in the forms required for the solver

function pooling_MIP_generation(initial_parameters::MIP_initial_parameters, generated_parameters::MIP_generated_parameters, Variables::Vector{pooling_variable}, bub, stoc_ppm::stochastic_pooling_problem_parameters )

    # Defining the model and the solver options
    original_problem = Model(() -> Gurobi.Optimizer(GRB_ENV))
    @suppress set_optimizer_attribute(original_problem, "NonConvex", initial_parameters.gurobi_parameters.NonConvex)
    @suppress set_optimizer_attribute(original_problem, "IntFeasTol",initial_parameters.gurobi_parameters.IntFeasTol)
    @suppress set_optimizer_attribute(original_problem, "FeasibilityTol", initial_parameters.gurobi_parameters.FeasibilityTol)
    @suppress set_optimizer_attribute(original_problem, "OptimalityTol", initial_parameters.gurobi_parameters.OptimalityTol)
    @suppress set_optimizer_attribute(original_problem, "Method", initial_parameters.gurobi_parameters.Method)
    #@suppress set_optimizer_attribute(original_problem, "OutputFlag", initial_parameters.gurobi_parameters.OutputFlag)
    @suppress set_optimizer_attribute(original_problem, "Threads", initial_parameters.gurobi_parameters.Threads)
    @suppress set_optimizer_attribute(original_problem, "TimeLimit", initial_parameters.solver_time_limit)
    #@suppress set_optimizer_attribute(original_problem, "Presolve", 0)

    # first stage decision variables 
    @variable(original_problem, x[ 1 : initial_parameters.num_first_stage_var, 1 : initial_parameters.num_scen ], binary = initial_parameters.bin_con_fs, integer = !initial_parameters.bin_con_fs)
    JuMP.unset_integer.(original_problem[:x][generated_parameters.x_cont_indexes, :])

    # second stage decision variables
    @variable(original_problem, y[ 1 : initial_parameters.num_second_stage_var, 1 : initial_parameters.num_scen ], Int)
    [JuMP.unset_integer.(original_problem[:y][generated_parameters.y_cont_indexes[:,j], j]) for j = 1:initial_parameters.num_scen]

    # slack variables
    @variable(original_problem, z[ 1 : initial_parameters.num_scen, 1 : initial_parameters.num_const ] >= 0 )

    # quadratic objective
    @objective(original_problem, Min,
        #map_coefficients_inplace!(a -> round(a, digits=1),
        sum( initial_parameters.scen_prob[s] *
            (
                sum( y[i, s] * generated_parameters.objective_Qs[s][i, j] * y[j, s] for i = 1 : initial_parameters.num_second_stage_var, j = 1 : initial_parameters.num_second_stage_var)
                + sum( x[i, s] * generated_parameters.objective_c[i] for i = 1 : initial_parameters.num_first_stage_var)
                + sum( y[i, s] * generated_parameters.objective_fs[s][i] for i = 1 : initial_parameters.num_second_stage_var)
                
            
            # slack variables 
            + initial_parameters.μ * sum(z[s,r] for r  = 1:initial_parameters.num_const ))
        for s in 1:initial_parameters.num_scen)
        )
        #)

    # quadratic constraints
    @constraint(original_problem, [s = 1 : initial_parameters.num_scen, i = 1 : initial_parameters.num_const ],
        sum( y[j, s] * generated_parameters.constraint_Qs[s][i][j,k] * y[k, s] for j = 1 : initial_parameters.num_second_stage_var, k = 1: initial_parameters.num_second_stage_var)
        + sum( x[j, s] * generated_parameters.constraint_fs[s][i][1, j] for j = 1 : initial_parameters.num_first_stage_var)
        + sum( y[j, s] * generated_parameters.constraint_fs[s][i][2, j] for j = 1 :  initial_parameters.num_second_stage_var)
        - generated_parameters.constraint_fs[s][i][3, 1] + z[s, i] == 0 )

    # box constraints for first stage variables
    @constraint(original_problem, [ s = 1 : initial_parameters.num_scen ],
        generated_parameters.x_boundaries[:, 1] .<= x[:, s] .<= generated_parameters.x_boundaries[:, 2])

    
    #---------------scenario-dependent constraints -----------------------

    # gathering the variables inidices that correpspond to quality at source nodes
    q_source = Variables[findall(x -> occursin("q[s", x.name) ,Variables)]

    #saving indices
    q_source_ind = Array{Int64}(undef,length(q_source))
    [q_source_ind[i] = q_source[i].index for i = 1:length(q_source)]

    # generate scenario-dependent source-node quality values
    Random.seed!(initial_parameters.random_seed) 
    q_source_value = Int.(round.(stoc_ppm.source_node_min_value .+ (stoc_ppm.source_node_max_value - stoc_ppm.source_node_min_value) .*  rand(length(q_source_ind), initial_parameters.num_scen), digits=0))

    # changing the boundary values to be the same for t=both left and right side for quality at source node 
    for s = 1: initial_parameters.num_scen
        for i = 1:length(q_source_ind)
            generated_parameters.y_boundaries[s][q_source_ind[i], 1] =  q_source_value[i]
            generated_parameters.y_boundaries[s][q_source_ind[i], 2] =  q_source_value[i]
        end
    end

    #---------------------------------------------------------------------------

    # box constraints for second stage variables 
    @constraint(original_problem, [s = 1 : initial_parameters.num_scen],
        generated_parameters.y_boundaries[s][:, 1] .<= y[:, s] .<= generated_parameters.y_boundaries[s][:, 2])

    # non-anticipativity conditions
    @constraint( original_problem, [s in 2 : initial_parameters.num_scen, i = 1 : initial_parameters.num_first_stage_var],
        x[i, s] - x[i, 1] == 0 )

    #---------------pooling problem related constraints -----------------------
    
    # flow balance
    number_of_pools = length(findall(x -> occursin("q[p", x.name), Variables))
    for i = 1:number_of_pools
        # gathering the variables inidices that correpspond to arcs to correspondent pool pi
        minus_arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) && occursin("p"*string(i)*"\")]", x.name) ,Variables)]
        # saving indices
        minus_arc_var_ind = Array{Int64}(undef,length(minus_arc_var))
        [minus_arc_var_ind[i] = minus_arc_var[i].index for i = 1:length(minus_arc_var)]

        # gathering the variables inidices that correpspond to arcs from correspondent pool pi
        plus_arc_var = Variables[findall(x -> occursin("x[Arc(\"p"*string(i), x.name), Variables)]
        # saving indices
        plus_arc_var_ind = Array{Int64}(undef,length(plus_arc_var))
        [plus_arc_var_ind[i] = plus_arc_var[i].index for i = 1:length(plus_arc_var)]
        
        @constraint(original_problem, [ s = 1 : initial_parameters.num_scen ],
            sum(y[minus_arc_var_ind, s]) - sum(y[plus_arc_var_ind, s]) == 0)   
    end


    # upper bounds of flows to target nodes 
    number_of_target_nodes = length(findall(x -> occursin("q[t", x.name), Variables))
    for i = 1:number_of_target_nodes
        # gathering the variables inidices that get to target node ti
        
        minus_arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) && occursin("t"*string(i)*"\")]", x.name) ,Variables)]
        # saving indices
        minus_arc_var_ind = Array{Int64}(undef,length(minus_arc_var))
        [minus_arc_var_ind[i] = minus_arc_var[i].index for i = 1:length(minus_arc_var)]
        
        @constraint(original_problem, [ s = 1 : initial_parameters.num_scen ],
            sum(y[minus_arc_var_ind, s]) <= bub[i])   
    end
    #---------------first-stage related constraints -----------------------
    # flow design
    # gathering the variables inidices that correpspond to all arcs
    arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) ,Variables)]
    
    #saving indices
    arc_var_ind = Array{Int64}(undef,length(arc_var))
    [arc_var_ind[i] = arc_var[i].index for i = 1:length(arc_var)]
    
    # making the flow constraint
    for i = 1:length(arc_var_ind)
        @constraint(original_problem, [ s = 1 : initial_parameters.num_scen ],
        y[arc_var_ind[i], s] <= stoc_ppm.flow_max_capacity * x[i,s] )
    end


    # pool deployment
    for i = 1:number_of_pools
        # gathering the variables inidices that correpspond to arcs to correspondent pool pi
        minus_arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) && occursin("p"*string(i)*"\")]", x.name) ,Variables)]
            
        # saving indices
        minus_arc_var_ind = Array{Int64}(undef,length(minus_arc_var))
        [minus_arc_var_ind[i] = minus_arc_var[i].index for i = 1:length(minus_arc_var)]

        @constraint(original_problem, [ s = 1 : initial_parameters.num_scen ],
            sum(y[minus_arc_var_ind, s]) <= stoc_ppm.pool_max_capacity * x[length(arc_var_ind) + i,s])       

    end

    #---------------scenario-dependent constraints -----------------------
    # lower bounds of flows to target nodes (demand)
    
    # generate scenario-dependent min demand values 
    Random.seed!(initial_parameters.random_seed) 
    min_demand = stoc_ppm.demand_min_value .+ (stoc_ppm.demand_max_value - stoc_ppm.demand_min_value) .*  rand(number_of_target_nodes, initial_parameters.num_scen)


    for i = 1:number_of_target_nodes
        # gathering the variables inidices that get to target node ti
        
        minus_arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) && occursin("t"*string(i)*"\")]", x.name) ,Variables)]
        # saving indices
        minus_arc_var_ind = Array{Int64}(undef,length(minus_arc_var))
        [minus_arc_var_ind[i] = minus_arc_var[i].index for i = 1:length(minus_arc_var)]
        
        @constraint(original_problem, [ s = 1 : initial_parameters.num_scen ],
            sum(y[minus_arc_var_ind, s]) >= min_demand[i,s])   
    end

return original_problem

end

function RNMDT_based_pooling_problem_generation(initial_parameters::MIP_initial_parameters, generated_parameters::MIP_generated_parameters, Variables::Vector{pooling_variable}, bub, stoc_ppm::stochastic_pooling_problem_parameters)

    # all the commments for the auxiliary constraints in this section refer to the paper
    # Enhancing the normalized multiparametric disaggregation
    # technique for mixed-integer quadratic programming

    # simplifying the notations
    precision_p = initial_parameters.RNMDT_precision_factor

    # defining RNMDT problem
    #RNMDT_problem = Model( optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol, "Method" => initial_parameters.gurobi_parameters.Method, "Threads" => initial_parameters.gurobi_parameters.Threads, "TimeLimit" => initial_parameters.solver_time_limit, "Presolve" => 0) )

    # defining RNMDT problem and solver options
    RNMDT_problem = Model(() -> Gurobi.Optimizer(GRB_ENV))
    @suppress set_optimizer_attribute(RNMDT_problem, "NonConvex", initial_parameters.gurobi_parameters.NonConvex)
    #@suppress set_optimizer_attribute(RNMDT_problem, "IntFeasTol",initial_parameters.gurobi_parameters.IntFeasTol)
    @suppress set_optimizer_attribute(RNMDT_problem, "FeasibilityTol", initial_parameters.gurobi_parameters.FeasibilityTol)
    @suppress set_optimizer_attribute(RNMDT_problem, "OptimalityTol", initial_parameters.gurobi_parameters.OptimalityTol)
    @suppress set_optimizer_attribute(RNMDT_problem, "Method", initial_parameters.gurobi_parameters.Method)
    #@suppress set_optimizer_attribute(RNMDT_problem, "OutputFlag", initial_parameters.gurobi_parameters.OutputFlag)
    @suppress set_optimizer_attribute(RNMDT_problem, "Threads", initial_parameters.gurobi_parameters.Threads)
    @suppress set_optimizer_attribute(RNMDT_problem, "TimeLimit", initial_parameters.solver_time_limit)
    #@suppress set_optimizer_attribute(RNMDT_problem, "Presolve", 0)

    # first stage decision variables
    @variable(RNMDT_problem, x[ 1 : initial_parameters.num_first_stage_var, 1 : initial_parameters.num_scen ], binary = initial_parameters.bin_con_fs, integer = !initial_parameters.bin_con_fs )
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
        sum( initial_parameters.scen_prob[s] *
            (
            sum(generated_parameters.objective_Qs[s][i, j] * w_RNMDT[i, j, s]
                    for i = 1 : initial_parameters.num_second_stage_var,
                        j = 1 : initial_parameters.num_second_stage_var)
                + sum( x[i, s] * generated_parameters.objective_c[i] for i = 1:initial_parameters.num_first_stage_var)
                + sum( y[j, s] * generated_parameters.objective_fs[s][j] for j = 1:initial_parameters.num_second_stage_var)
            # slack variables
            + initial_parameters.μ * sum(z[s,r] for r  = 1:initial_parameters.num_const )          
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
        - generated_parameters.constraint_fs[s][r][3, 1] + z[s,r] == 0
    ) # 27

    ## auxiliary RNMDT constraints
    ss_nf  = Array(1 : initial_parameters.num_second_stage_var)
    #---------------scenario-dependent constraints -----------------------

    # gathering the variables inidices that correpspond to quality at source nodes
    q_source = Variables[findall(x -> occursin("q[s", x.name) ,Variables)]

    #saving indices
    q_source_ind = Array{Int64}(undef,length(q_source))
    [q_source_ind[i] = q_source[i].index for i = 1:length(q_source)]

    #@show initial_parameters.num_second_stage_var

    # generate scenario-dependent source-node quality values
    Random.seed!(initial_parameters.random_seed) 
    q_source_value = Int.(round.(stoc_ppm.source_node_min_value .+ (stoc_ppm.source_node_max_value - stoc_ppm.source_node_min_value) .*  rand(length(q_source_ind), initial_parameters.num_scen), digits=0))

    # changing the boundary values to be the same for t=both left and right side for qualities at source nodes 
    for s = 1: initial_parameters.num_scen
        for i = 1:length(q_source_ind)
            fix(y[q_source_ind[i], s], q_source_value[i]; force = true)
            for j = 1:length(q_source_ind)
                fix( w_RNMDT[q_source_ind[i], q_source_ind[j], s], q_source_value[i]*q_source_value[j]; force = true)
            end
        end
    end

    # changing the boundary values to be the same for t=both left and right side for quality at source node 
    for s = 1: initial_parameters.num_scen
           for i = 1:length(q_source_ind)
               #generated_parameters.y_boundaries[s][q_source_ind[i], 1] =  q_source_value[i]
               generated_parameters.y_boundaries[s][q_source_ind[i], 1] = 0
               generated_parameters.y_boundaries[s][q_source_ind[i], 2] =  q_source_value[i]
           end
    end

    deleteat!(ss_nf, q_source_ind)
    #@show ss_nf
    
    #---------------------------------------------------------------------------
   
    @constraint( RNMDT_problem, [ j  in 1 : initial_parameters.num_second_stage_var, s  = 1 : initial_parameters.num_scen],
        y[j, s] == (generated_parameters.y_boundaries[s][j, 2] - generated_parameters.y_boundaries[s][j, 1])  *
        ( sum( 2.0^l * z_RNMDT[j, l, s] for l = precision_p[j, s] : -1) + delta_y_RNMDT[j, s] ) ) # 28

    @constraint( RNMDT_problem, [ i in 1 : initial_parameters.num_second_stage_var, j in 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen ],
        w_RNMDT[i, j, s] == (generated_parameters.y_boundaries[s][j, 2] - generated_parameters.y_boundaries[s][j, 1]) *
        ( sum( 2.0^l * y_heat_RNMDT[i, j, l, s] for l = precision_p[j, s] : -1) + delta_w_RNMDT[i, j, s] ) ) # 29

    @constraint( RNMDT_problem, [ j  in 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen ],
        0 <= delta_y_RNMDT[j, s] <= 2.0^precision_p[j, s] ) # 30

    @constraint( RNMDT_problem, [ i in 1 : initial_parameters.num_second_stage_var, j in 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen],
        2.0^precision_p[i, s] * ( y[i, s] - generated_parameters.y_boundaries[s][i, 2] ) + generated_parameters.y_boundaries[s][i, 2] * delta_y_RNMDT[j, s] <= delta_w_RNMDT[i, j, s]) # 31

    @constraint( RNMDT_problem, [ i in 1 : initial_parameters.num_second_stage_var, j in 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen],
        delta_w_RNMDT[i, j, s] <= 2.0^precision_p[i, s] * ( y[i, s] - generated_parameters.y_boundaries[s][i, 1] ) + generated_parameters.y_boundaries[s][i, 1] * delta_y_RNMDT[j, s]) # 31

    @constraint( RNMDT_problem, [ i in 1 : initial_parameters.num_second_stage_var, j in 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen ],
        generated_parameters.y_boundaries[s][i, 1] * delta_y_RNMDT[j, s] <= delta_w_RNMDT[i, j, s]) #32

    @constraint( RNMDT_problem, [ i in 1 : initial_parameters.num_second_stage_var, j in 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen],
        delta_w_RNMDT[i, j, s] <= generated_parameters.y_boundaries[s][i, 2] * delta_y_RNMDT[j, s]) # 32

    @constraint( RNMDT_problem, [ i in 1 : initial_parameters.num_second_stage_var, j in 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen, l = precision_p[j, s] : -1 ],
        generated_parameters.y_boundaries[s][i, 1] * z_RNMDT[j, l, s]  <= y_heat_RNMDT[i, j, l, s]) # 33

    @constraint( RNMDT_problem, [ i in 1 : initial_parameters.num_second_stage_var, j in 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen, l = precision_p[j, s] : -1 ],
        y_heat_RNMDT[i, j, l, s] <= generated_parameters.y_boundaries[s][i, 2] * z_RNMDT[j, l, s]) # 33

    @constraint( RNMDT_problem, [ i in 1 : initial_parameters.num_second_stage_var, j in 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen, l = precision_p[j, s] : -1 ],
        generated_parameters.y_boundaries[s][i, 1] * (1 - z_RNMDT[j, l, s]) <= y[i,s] - y_heat_RNMDT[i, j, l, s] ) # 34

    @constraint( RNMDT_problem, [ i in 1 : initial_parameters.num_second_stage_var, j in 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen, l = precision_p[j, s] : -1 ],
        y[i,s] - y_heat_RNMDT[i, j, l, s] <= generated_parameters.y_boundaries[s][i, 2] * (1 - z_RNMDT[j, l, s]) ) # 34

    ##
    # box constraints for first stage variables
    @constraint(RNMDT_problem, [ s = 1 : initial_parameters.num_scen ],
        generated_parameters.x_boundaries[:, 1] .<= x[:, s] .<= generated_parameters.x_boundaries[:, 2])
    

    # box constraints for second stage variables
    @constraint(RNMDT_problem, [i in 1 : initial_parameters.num_second_stage_var, s = 1 : initial_parameters.num_scen],
        generated_parameters.y_boundaries[s][i, 1] .<= y[i, s] .<= generated_parameters.y_boundaries[s][i, 2])

    # non-anticipativity conditions
    @constraint(RNMDT_problem, [s in 2 : initial_parameters.num_scen, i = 1 : initial_parameters.num_first_stage_var],
        x[i, s] - x[i, 1] == 0 )

    #---------------pooling problem related constraints -----------------------
    
    # flow balance
    number_of_pools = length(findall(x -> occursin("q[p", x.name), Variables))
    for i = 1:number_of_pools
        # gathering the variables inidices that correpspond to arcs to correspondent pool pi
        minus_arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) && occursin("p"*string(i)*"\")]", x.name) ,Variables)]
        # saving indices
        minus_arc_var_ind = Array{Int64}(undef,length(minus_arc_var))
        [minus_arc_var_ind[i] = minus_arc_var[i].index for i = 1:length(minus_arc_var)]

        # gathering the variables inidices that correpspond to arcs from correspondent pool pi
        plus_arc_var = Variables[findall(x -> occursin("x[Arc(\"p"*string(i), x.name), Variables)]
        # saving indices
        plus_arc_var_ind = Array{Int64}(undef,length(plus_arc_var))
        [plus_arc_var_ind[i] = plus_arc_var[i].index for i = 1:length(plus_arc_var)]
        
        @constraint(RNMDT_problem, [ s = 1 : initial_parameters.num_scen ],
            sum(y[minus_arc_var_ind, s]) - sum(y[plus_arc_var_ind, s]) == 0)   
    end


    # upper bounds of flows to target nodes 
    number_of_target_nodes = length(findall(x -> occursin("q[t", x.name), Variables))
    for i = 1:number_of_target_nodes
        # gathering the variables inidices that get to target node ti
        
        minus_arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) && occursin("t"*string(i)*"\")]", x.name) ,Variables)]
        # saving indices
        minus_arc_var_ind = Array{Int64}(undef,length(minus_arc_var))
        [minus_arc_var_ind[i] = minus_arc_var[i].index for i = 1:length(minus_arc_var)]
        
        @constraint(RNMDT_problem, [ s = 1 : initial_parameters.num_scen ],
            sum(y[minus_arc_var_ind, s]) <= bub[i])   
    end
    #---------------first-stage related constraints -----------------------
    # flow design
    # gathering the variables inidices that correpspond to all arcs
    arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) ,Variables)]
    
    #saving indices
    arc_var_ind = Array{Int64}(undef,length(arc_var))
    [arc_var_ind[i] = arc_var[i].index for i = 1:length(arc_var)]
    
    # making the flow constraint
    for i = 1:length(arc_var_ind)
        @constraint(RNMDT_problem, [ s = 1 : initial_parameters.num_scen ],
            y[arc_var_ind[i], s] <= stoc_ppm.flow_max_capacity * x[i,s] )
    end


    # pool deployment
    for i = 1:number_of_pools
        # gathering the variables inidices that correpspond to arcs to correspondent pool pi
        minus_arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) && occursin("p"*string(i)*"\")]", x.name) ,Variables)]
            
        # saving indices
        minus_arc_var_ind = Array{Int64}(undef,length(minus_arc_var))
        [minus_arc_var_ind[i] = minus_arc_var[i].index for i = 1:length(minus_arc_var)]

        @constraint(RNMDT_problem, [ s = 1 : initial_parameters.num_scen ],
            sum(y[minus_arc_var_ind, s]) <= stoc_ppm.pool_max_capacity * x[length(arc_var_ind) + i,s])       

    end

    #---------------scenario-dependent constraints -----------------------
    # lower bounds of flows to target nodes (demand)
    
    # generate scenario-dependent min demand values 
    Random.seed!(initial_parameters.random_seed) 
    min_demand = stoc_ppm.demand_min_value .+ (stoc_ppm.demand_max_value - stoc_ppm.demand_min_value) .*  rand(number_of_target_nodes, initial_parameters.num_scen)


    for i = 1:number_of_target_nodes
        # gathering the variables inidices that get to target node ti
        
        minus_arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) && occursin("t"*string(i)*"\")]", x.name) ,Variables)]
        # saving indices
        minus_arc_var_ind = Array{Int64}(undef,length(minus_arc_var))
        [minus_arc_var_ind[i] = minus_arc_var[i].index for i = 1:length(minus_arc_var)]
        
        @constraint(RNMDT_problem, [ s = 1 : initial_parameters.num_scen ],
            sum(y[minus_arc_var_ind, s]) >= min_demand[i,s])   
    end

    return RNMDT_problem
end


function RNMDT_based_augmented_lagrangian_relaxation_pooling_problem_generation(initial_parameters::MIP_initial_parameters, generated_parameters::MIP_generated_parameters, Variables::Vector{pooling_variable}, bub, stoc_ppm::stochastic_pooling_problem_parameters)

    # simplifying the notations
    precision_p = initial_parameters.RNMDT_precision_factor

    # generate the starting values for the lagrangian multipliers
    vector_of_lambda_lagrangian = Array{Array{Float64}}(undef, initial_parameters.num_scen)
    [ vector_of_lambda_lagrangian[i] = zeros(1,initial_parameters.num_first_stage_var)
            for i = 1 : initial_parameters.num_scen ]

    # creating the array of subproblems
    vector_of_subproblems = Array{Any}(undef, initial_parameters.num_scen)

    # formulating the subproblems based on the scenarios
    for s = 1 : initial_parameters.num_scen

        #vector_of_subproblems[s] = Model(optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => initial_parameters.gurobi_parameters.NonConvex, "IntFeasTol" =>  initial_parameters.gurobi_parameters.IntFeasTol, "FeasibilityTol" =>  initial_parameters.gurobi_parameters.FeasibilityTol, "OptimalityTol" =>  initial_parameters.gurobi_parameters.OptimalityTol, "Method" => initial_parameters.gurobi_parameters.Method, "OutputFlag" => initial_parameters.gurobi_parameters.OutputFlag,  "Threads" => #initial_parameters.gurobi_parameters.Threads, "NumericFocus" => initial_parameters.gurobi_parameters.NumericFocus))

        # defining subproblems and solver options
        vector_of_subproblems[s] = Model(() -> Gurobi.Optimizer(GRB_ENV))
        @suppress set_optimizer_attribute(vector_of_subproblems[s], "NonConvex", initial_parameters.gurobi_parameters.NonConvex)
        @suppress set_optimizer_attribute(vector_of_subproblems[s], "IntFeasTol",initial_parameters.gurobi_parameters.IntFeasTol)
        @suppress set_optimizer_attribute(vector_of_subproblems[s], "FeasibilityTol", initial_parameters.gurobi_parameters.FeasibilityTol)
        @suppress set_optimizer_attribute(vector_of_subproblems[s], "OptimalityTol", initial_parameters.gurobi_parameters.OptimalityTol)
        #@suppress set_optimizer_attribute(vector_of_subproblems[s], "Method", initial_parameters.gurobi_parameters.Method)
        #@suppress set_optimizer_attribute(vector_of_subproblems[s], "OutputFlag", initial_parameters.gurobi_parameters.OutputFlag)
        @suppress set_optimizer_attribute(vector_of_subproblems[s], "Threads", initial_parameters.gurobi_parameters.Threads)
        #@suppress set_optimizer_attribute(vector_of_subproblems[s], "NumericFocus", initial_parameters.gurobi_parameters.NumericFocus)
        #@suppress set_optimizer_attribute(vector_of_subproblems[s], "Presolve", 0)

        # first stage decision variables
        @variable(vector_of_subproblems[s], x[1 : initial_parameters.num_first_stage_var], binary = initial_parameters.bin_con_fs, integer = !initial_parameters.bin_con_fs )
        JuMP.unset_integer.(vector_of_subproblems[s][:x][generated_parameters.x_cont_indexes])

        # non-anticipativity conditions variable (augmented lagrangian)
        @variable(vector_of_subproblems[s], al_z[1 : initial_parameters.num_first_stage_var]  )
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
            (initial_parameters.scen_prob[s] *
            ( sum(generated_parameters.objective_Qs[s][i, j] * w_RNMDT[i, j]
                for i = 1 : initial_parameters.num_second_stage_var,
                    j = 1 : initial_parameters.num_second_stage_var)
            + sum( x[i] * generated_parameters.objective_c[i]  for i = 1:initial_parameters.num_first_stage_var)
            + sum( y[j] * generated_parameters.objective_fs[s][j]  for j = 1:initial_parameters.num_second_stage_var)
            #slack variables 
            + initial_parameters.μ * sum(z[r] for r  = 1:initial_parameters.num_const )
            )
            + sum( vector_of_lambda_lagrangian[s] .* (x .- al_z) )

            + sum( initial_parameters.al_penalty_parameter/2 .* (x .- al_z) .* (x .- al_z) )
            )

        )

        #quadratic constraint
        @constraint( vector_of_subproblems[s], [ r = 1 : initial_parameters.num_const ],
            sum(generated_parameters.constraint_Qs[s][r][i, j] * w_RNMDT[i, j]
                for i = 1 : initial_parameters.num_second_stage_var,
                j = 1 : initial_parameters.num_second_stage_var)
            + sum( x[i] * generated_parameters.constraint_fs[s][r][1, i] for i = 1:initial_parameters.num_first_stage_var)
            + sum( y[j] * generated_parameters.constraint_fs[s][r][2, j] for j = 1:initial_parameters.num_second_stage_var )
            - generated_parameters.constraint_fs[s][r][3, 1] + z[r] == 0
        ) # 27

        # box constraints for first stage variables
        @constraint( vector_of_subproblems[s],
            generated_parameters.x_boundaries[:, 1] .<= x .<= generated_parameters.x_boundaries[:, 2])

        #---------------scenario-dependent constraints -----------------------

        # gathering the variables inidices that correpspond to quality at source nodes
        q_source = Variables[findall(x -> occursin("q[s", x.name) ,Variables)]

        #saving indices
        q_source_ind = Array{Int64}(undef,length(q_source))
        [q_source_ind[i] = q_source[i].index for i = 1:length(q_source)]

        # generate scenario-dependent source-node quality values
        Random.seed!(initial_parameters.random_seed) 
        q_source_value = Int.(round.(stoc_ppm.source_node_min_value .+ (stoc_ppm.source_node_max_value - stoc_ppm.source_node_min_value) .*  rand(length(q_source_ind), initial_parameters.num_scen), digits=0))

        # changing the boundary values to be the same for t=both left and right side for quality at source node 
        for i = 1:length(q_source_ind)
            generated_parameters.y_boundaries[s][q_source_ind[i], 1] =  q_source_value[i]
            generated_parameters.y_boundaries[s][q_source_ind[i], 2] =  q_source_value[i]
        end

        #---------------------------------------------------------------------------

        # box constraints for second stage variables
        @constraint( vector_of_subproblems[s],
            generated_parameters.y_boundaries[s][:, 1] .<= y .<= generated_parameters.y_boundaries[s][:, 2])

    ## auxiliary RNMDT constraints
        

        #---------------scenario-dependent constraints -----------------------
        ss_nf  = Array(1 : initial_parameters.num_second_stage_var)

        # gathering the variables inidices that correpspond to quality at source nodes
        q_source = Variables[findall(x -> occursin("q[s", x.name) ,Variables)]

        #saving indices
        q_source_ind = Array{Int64}(undef,length(q_source))
        [q_source_ind[i] = q_source[i].index for i = 1:length(q_source)]

        #@show initial_parameters.num_second_stage_var

        # generate scenario-dependent source-node quality values
        Random.seed!(initial_parameters.random_seed) 
        q_source_value = Int.(round.(stoc_ppm.source_node_min_value .+ (stoc_ppm.source_node_max_value - stoc_ppm.source_node_min_value) .*  rand(length(q_source_ind), initial_parameters.num_scen), digits=0))

        # changing the boundary values to be the same for t=both left and right side for qualities at source nodes 
        for i = 1:length(q_source_ind)
            fix(y[q_source_ind[i]], q_source_value[i]; force = true)
            for j = 1:length(q_source_ind)
                fix( w_RNMDT[q_source_ind[i], q_source_ind[j]], q_source_value[i]*q_source_value[j]; force = true)
            end
        end

        deleteat!(ss_nf, q_source_ind)
        #@show ss_nf
        #---------------------------------------------------------------------------

        @constraint( vector_of_subproblems[s], [ j in  ss_nf],
            y[j] == (generated_parameters.y_boundaries[s][j, 2] - generated_parameters.y_boundaries[s][j, 1])  *
            ( sum( 2.0^l * z_RNMDT[j, l] for l = precision_p[j, s] : -1) + delta_y_RNMDT[j] ) ) # 28

        @constraint( vector_of_subproblems[s], [ i  in  ss_nf, j in  ss_nf ],
            w_RNMDT[i, j] == (generated_parameters.y_boundaries[s][j, 2] - generated_parameters.y_boundaries[s][j, 1]) *
            ( sum( 2.0^l * y_heat_RNMDT[i, j, l] for l = precision_p[j, s] : -1) + delta_w_RNMDT[i, j] ) ) # 29

        @constraint( vector_of_subproblems[s], [ j  in  ss_nf ],
            0 <= delta_y_RNMDT[j] <= 2.0^precision_p[j, s] ) # 30

        @constraint( vector_of_subproblems[s], [ i in  ss_nf, j in  ss_nf],
            2.0^precision_p[i, s] * ( y[i] - generated_parameters.y_boundaries[s][i, 2] ) + generated_parameters.y_boundaries[s][i, 2] * delta_y_RNMDT[j] <= delta_w_RNMDT[i, j]) # 31

        @constraint( vector_of_subproblems[s], [ i in  ss_nf, j in  ss_nf],
            delta_w_RNMDT[i, j] <= 2.0^precision_p[i, s] * ( y[i] - generated_parameters.y_boundaries[s][i, 1] ) + generated_parameters.y_boundaries[s][i, 1] * delta_y_RNMDT[j]) # 31

        @constraint( vector_of_subproblems[s], [ i in  ss_nf, j in  ss_nf],
            generated_parameters.y_boundaries[s][i, 1] * delta_y_RNMDT[j] <= delta_w_RNMDT[i, j]) #32

        @constraint( vector_of_subproblems[s], [ i in  ss_nf, j in  ss_nf],
            delta_w_RNMDT[i, j] <= generated_parameters.y_boundaries[s][i, 2] * delta_y_RNMDT[j]) # 32

        @constraint( vector_of_subproblems[s], [ i in  ss_nf, j in  ss_nf, l = precision_p[j, s] : -1],
            generated_parameters.y_boundaries[s][i, 1] * z_RNMDT[j, l]  <= y_heat_RNMDT[i, j, l]) # 33

        @constraint( vector_of_subproblems[s], [ i in  ss_nf, j in  ss_nf, l = precision_p[j, s] : -1],
            y_heat_RNMDT[i, j, l] <= generated_parameters.y_boundaries[s][i, 2] * z_RNMDT[j, l]) # 33

        @constraint( vector_of_subproblems[s], [ i in  ss_nf, j in  ss_nf, l = precision_p[j, s] : -1],
            generated_parameters.y_boundaries[s][i, 1] * (1 - z_RNMDT[j, l]) <= y[i] - y_heat_RNMDT[i, j, l] ) # 34

        @constraint( vector_of_subproblems[s], [ i in  ss_nf, j in  ss_nf, l = precision_p[j, s] : -1],
            y[i] - y_heat_RNMDT[i, j, l] <= generated_parameters.y_boundaries[s][i, 2] * (1 - z_RNMDT[j, l]) ) # 34
    ##  #---------------pooling problem related constraints -----------------------
    
        # flow balance
        number_of_pools = length(findall(x -> occursin("q[p", x.name), Variables))
        for i = 1:number_of_pools
            # gathering the variables inidices that correpspond to arcs to correspondent pool pi
            minus_arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) && occursin("p"*string(i)*"\")]", x.name) ,Variables)]
            # saving indices
            minus_arc_var_ind = Array{Int64}(undef,length(minus_arc_var))
            [minus_arc_var_ind[i] = minus_arc_var[i].index for i = 1:length(minus_arc_var)]

            # gathering the variables inidices that correpspond to arcs from correspondent pool pi
            plus_arc_var = Variables[findall(x -> occursin("x[Arc(\"p"*string(i), x.name), Variables)]
            # saving indices
            plus_arc_var_ind = Array{Int64}(undef,length(plus_arc_var))
            [plus_arc_var_ind[i] = plus_arc_var[i].index for i = 1:length(plus_arc_var)]
        
            @constraint(vector_of_subproblems[s],
                sum(y[minus_arc_var_ind]) - sum(y[plus_arc_var_ind]) == 0)   
        end


        # upper bounds of flows to target nodes 
        number_of_target_nodes = length(findall(x -> occursin("q[t", x.name), Variables))
        for i = 1:number_of_target_nodes
            # gathering the variables inidices that get to target node ti
        
            minus_arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) && occursin("t"*string(i)*"\")]", x.name) ,Variables)]
            # saving indices
            minus_arc_var_ind = Array{Int64}(undef,length(minus_arc_var))
            [minus_arc_var_ind[i] = minus_arc_var[i].index for i = 1:length(minus_arc_var)]
        
            @constraint(vector_of_subproblems[s],
                sum(y[minus_arc_var_ind]) <= bub[i])   
        end
        #---------------first-stage related constraints -----------------------
        # flow design
        # gathering the variables inidices that correpspond to all arcs
        arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) ,Variables)]
    
        #saving indices
        arc_var_ind = Array{Int64}(undef,length(arc_var))
        [arc_var_ind[i] = arc_var[i].index for i = 1:length(arc_var)]
    
        # making the flow constraint
        for i = 1:length(arc_var_ind)
            @constraint(vector_of_subproblems[s],
                y[arc_var_ind[i]] <= stoc_ppm.flow_max_capacity * x[i] )
        end


        # pool deployment
        for i = 1:number_of_pools
            # gathering the variables inidices that correpspond to arcs to correspondent pool pi
            minus_arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) && occursin("p"*string(i)*"\")]", x.name) ,Variables)]
            
            # saving indices
            minus_arc_var_ind = Array{Int64}(undef,length(minus_arc_var))
            [minus_arc_var_ind[i] = minus_arc_var[i].index for i = 1:length(minus_arc_var)]

            @constraint(vector_of_subproblems[s],
                sum(y[minus_arc_var_ind]) <= stoc_ppm.pool_max_capacity * x[length(arc_var_ind) + i])       

        end

        #---------------scenario-dependent constraints -----------------------
        # lower bounds of flows to target nodes (demand)
    
        # generate scenario-dependent min demand values 
        Random.seed!(initial_parameters.random_seed) 
        min_demand = stoc_ppm.demand_min_value .+ (stoc_ppm.demand_max_value - stoc_ppm.demand_min_value) .*  rand(number_of_target_nodes, initial_parameters.num_scen)


        for i = 1:number_of_target_nodes
            # gathering the variables inidices that get to target node ti
        
            minus_arc_var = Variables[findall(x -> occursin("x[Arc(", x.name) && occursin("t"*string(i)*"\")]", x.name) ,Variables)]
            # saving indices
            minus_arc_var_ind = Array{Int64}(undef,length(minus_arc_var))
            [minus_arc_var_ind[i] = minus_arc_var[i].index for i = 1:length(minus_arc_var)]
        
            @constraint(vector_of_subproblems[s],
                sum(y[minus_arc_var_ind]) >= min_demand[i,s])   
        end

    end

    return vector_of_subproblems

end