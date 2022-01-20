src_link  =  "/scratch/work/belyakn1/BnB_p_lagrangian/src/"
#src_link  =  "/Users/nikitabelyak/Dropbox (Aalto)/branch-and-bound-caroe-and-schultz/src/"

cd(src_link)
using Pkg
Pkg.activate(".")
Pkg.update()
#Pkg.resolve()
#ENV["GUROBI_HOME"]=  "/Library/gurobi910/mac64/"

ENV["GUROBI_HOME"] = "/share/apps/anaconda-ci/fgci-centos7-generic/software/anaconda/2020-01-tf2/5a34a04a"
ENV["GRB_LICENSE_FILE"]="/scratch/work/belyakn1/BnB_p_lagrangian/gurobi.lic"
Pkg.build("Gurobi")
Pkg.instantiate()

include(src_link*"triton_initialization.jl")

const GRB_ENV = Gurobi.Env()


# check whether the folder for the "today" experiments exists
# and if not create one
if !isdir(chop(src_link, tail = 4) * "experiments_" * string(Dates.today()))
    mkdir(chop(src_link, tail = 4) * "experiments_" * string(Dates.today()))
end

output_link = chop(src_link, tail = 4) * "experiments_" * string(Dates.today()) * "/"

## Constructing the experiments
# the structure that will collect the experiments results
#output_df = DataFrame( num_of_scen = Int[], num_fs_var = Int[], num_ss_var = Int[], num_const = Int[], p_RNMDT = Int[], primal_f = Float64[], primal_x = String[], primal_gap = Float64[], RNMDT_UB = Float64[], RNMDT_x = String[], RNMDT_time = Float64[], RNMDT_wy_gap = Float64[], BnB_UB = Float64[], BnB_LB = Float64[], BnB_x = String[], BnB_time = Float64[], BnB_wy_gap = Float64[], BnB_nodes_explored = Int[] )

## Defining the initial parameters for the experiments

# the number of scenarios
n_scenarios = [10, 15, 20]
# the number of the first stage variables / per scenario
fs_var = [10, 15, 20]
# the number of the second stage variabes / per scenario
ss_var = [15, 20, 25]
# the number of the constraints / per scenario
const_num = [15, 20, 25]
# Methods to be used [full scale, RNMDT, BnB]: 1 = use, 0 = don't use
experiments_methods = [1, 1, 1]
# The maximum time limit for the set of instances
g_time_limit = 3600
# precision factor values to be considered for RNMDT
p_value = [-2 -4 -6 -8]

## Constructing the experiments
# the structure that will collect the experiments results

function experiments_function(n_scenarios, fs_var, ss_var, const_num, experiments_methods, g_time_limit, p_values, output_link )

    # if we solve the instances using all the three methods
    if sum(experiments_methods) == length(experiments_methods)

        # forming corresponding data frame
        output_df = DataFrame( num_of_scen = Int[], num_fs_var = Int[], num_ss_var = Int[], num_const = Int[], p_RNMDT = Int[], primal_f = Float64[], primal_x = String[], primal_gap = Float64[], RNMDT_UB = Float64[], RNMDT_x = String[], RNMDT_time = Float64[], RNMDT_wy_gap = Float64[], BnB_UB = Float64[], BnB_LB = Float64[], BnB_x = String[], BnB_time = Float64[], BnB_wy_gap = Float64[], BnB_nodes_explored = Int[] )

        for s in n_scenarios
            for i_fs_var in fs_var
                for i_ss_var in ss_var
                    for i_const_num in const_num
                        
                        # defining auxiliary variable to keep us track at which index we are at the p_values array 
                        p_index = 1
                        
                        # defining starting value of the precision factor
                        p = p_values[p_index]

                        # defining general intial parameters for the optimisation problem
                        initial_parameters = initialisation(s,i_fs_var,i_ss_var,i_const_num, p)

                        # generating the structure containing the constraints and objective related parameters
                        generated_parameters = parameters_generation(initial_parameters)

                        # soving the full scale problem using global solver
                        primal_problem = MIP_generation(initial_parameters, generated_parameters)
                        @suppress optimize!(primal_problem)
                        primal_problem_f = objective_value(primal_problem)
                        primal_problem_x = string(value.(primal_problem[:x][:,1]))
                        primal_problem_optimality_gap = MOI.get(primal_problem, MOI.RelativeGap())

                        # defining the variables corresponing to the total time spent by the global solver
                        # and developed method for solving the RNMDT isntances with different values of the
                        # precision factor p
                        bnb_g_time = 0
                        rnmdt_g_time = 0

                        # forming an array to collect the dual objective values resutling from iterations of FW-Ph
                        # applied to the instances with different values of the precision factor p
                        FW_PH_dual_objective_output = []

                        # forming an array to collect the info whether we solved the RNMDT instance with the global solver or not
                        RNMDT_global_solved = []

                        # collecting the objective values obtained by solving RNMDT instance with the global solver
                        RNMDT_global_objective = []

                        # the variable to keep the info on the max number of iterations used by FW-PH
                        # for the "fair" plot
                        FW_PH_max_iter_used = 0


                        # here we only consider BnB-related global time since it is our primal focus
                        # and it is generally assumed to be faster => allowing for solving higer of number of instances
                        # while precision factor is decreasing
                        while (bnb_g_time < g_time_limit) && (p_index <= length(p_values))
                            
                            if p_index > 1

                                p = p_values[p_index]

                                # defining general intial parameters for the optimisation problem
                                initial_parameters = initialisation(s,i_fs_var,i_ss_var,i_const_num, p)

                                # generating the structure containing the constraints and objective related parameters
                                generated_parameters = parameters_generation(initial_parameters)
                            
                            end
                            
                            RNMDT_relaxation = RNMDT_based_problem_generation(initial_parameters, generated_parameters)

                            # if we did not exceed yet the overall time allowed for solving instances
                            # with different values of the precision factor p
                            if rnmdt_g_time < g_time_limit
                                # solving the RNMDT instance using the global solver
                                rndmt_init_time = time()
                                @suppress optimize!(RNMDT_relaxation)
                                rndmt_final_time = time() - rndmt_init_time
                                # updating the total time spent by the global solver to solve RNMDT instances
                                rnmdt_g_time += rndmt_final_time
                                #rndmt_happened = true
                                push!(RNMDT_global_solved, true)
                                push!(RNMDT_global_objective, objective_value(RNMDT_relaxation))
                            else
                                #rndmt_happened = false
                                push!(RNMDT_global_solved, false)
                            end # if

                            # if we need to plot the FW iterations for the root node
                            if initial_parameters.PH_SDM_parameters.PH_plot

                                # solving the RNMDT instance using the BnB+FW-PH method
                                bnb_p_init_time = time()
                                bnb_output, bnb_FW_iterations_output = bnb_solve(initial_parameters, non_ant_tol, tol_bb, integrality_tolerance)
                                bnb_p_final_time = time() - bnb_p_init_time

                                # collecting the data for plot
                                push!(FW_PH_dual_objective_output,bnb_FW_iterations_output[1])

                                # updating the max_number of the iteratiosn used by FW-PH for the different values of the precision factor p
                                length(bnb_FW_iterations_output[1])>FW_PH_max_iter_used ? FW_PH_max_iter_used = length(bnb_FW_iterations_output[1]) : true

                                # froming the dual_feasibility related output to be suitable for DataFrame
                                dual_feasibility_output = []
                                for i = 1:length(bnb_FW_iterations_output[2])
                                    push!(dual_feasibility_output, string(bnb_FW_iterations_output[2][i]))
                                end
                                # Printing the dual_deasibility_condtion and primal_dual_residual resulting from iterations of FW-PH applied to the root node
                                FW_PH_output = DataFrame( iteration = 1:length(bnb_FW_iterations_output[2]), dual_deasibility_cond = dual_feasibility_output, primal_dual_residual = bnb_FW_iterations_output[3])
                                XLSX.writetable(output_link*"$(s)_$(i_fs_var)_$(i_ss_var)_$(i_const_num)_$(p)_" * string(Dates.now()) * "Dual feas cond and primal dual res " * ".xlsx", FW_PH_output)

                                # printing out the info on the identical etries in the feasibility set V
                                io1 = open(output_link*"$(s)_$(i_fs_var)_$(i_ss_var)_$(i_const_num)_$(p)_" * string(Dates.now())* " feasibility set V identical entries count" * ".txt", "w")
                                for s = 1:length(bnb_FW_iterations_output[5])
                                        println(io1, "\n SCENARIO = $s : $(bnb_FW_iterations_output[5][s])")
                                end
                                close(io1)

                                # Printing the fesibility set V resulting from each iteration of FW-PH applied to the root node
                                io = open(output_link*"$(s)_$(i_fs_var)_$(i_ss_var)_$(i_const_num)_$(p)_" * string(Dates.now())* " feasibility set V " * ".txt", "w")
                                    # printing only for the first scenario to ease representation
                                    for i = 1:length(bnb_FW_iterations_output[4])
                                        println(io, "\n ITERATION = $i : $(length(bnb_FW_iterations_output[4][i][1]))")
                                    end
                                    println(io, "--------------------------------------------------------")
                                    for i = 1:length(bnb_FW_iterations_output[4])
                                        println(io, "\n ITERATION = $i \n")
                                        println(io, bnb_FW_iterations_output[4][i][1])
                                        println(io,"\n")
                                    end
                                close(io)
                            else
                                # solving the RNMDT instance using the BnB+FW-PH method
                                bnb_p_init_time = time()
                                bnb_output = bnb_solve(initial_parameters, non_ant_tol, tol_bb, integrality_tolerance)
                                bnb_p_final_time = time() - bnb_p_init_time
                            end



                            bnb_g_time += bnb_p_final_time

                            # adding the line to the data base with the RNMDt-related information correspondent to whether
                            # we have solved the RNMDT instance for the current value of precision factor p with the
                            # global solver or not
                            if RNMDT_global_solved[end]
                                push!(output_df, (s, i_fs_var, i_ss_var, i_const_num , p, primal_problem_f, primal_problem_x, primal_problem_optimality_gap, objective_value(RNMDT_relaxation), string(value.(RNMDT_relaxation[:x][:,1])), rndmt_final_time, RNMDT_gap_computation( value.(RNMDT_relaxation[:y]), value.(RNMDT_relaxation[:w_RNMDT])), bnb_output.UBDg, bnb_output.LBDg, string(bnb_output.soln_val), bnb_p_final_time, bnb_output.RNMDT_gap_wy, bnb_output.nodes_used))
                            else
                                push!(output_df, (s, i_fs_var, i_ss_var, i_const_num, p, primal_problem_f, primal_problem_x, primal_problem_optimality_gap, 0.0, "NaN", 0.0, 0.0, bnb_output.UBDg, bnb_output.LBDg, string(bnb_output.soln_val), bnb_p_final_time,  bnb_output.RNMDT_gap_wy, bnb_output.nodes_used ))
                            end # if

                            # moving forward among the values of the precision factor 
                            p_index += 1

                            # saving the updated database after solving each instance (in case there are some technical issues stopping the code)
                            XLSX.writetable(output_link*"experiments"*string(Dates.now())*".xlsx", output_df)
                        end # while

                        # ploting the FW-PH dual objective related output
                        for i = 1:length(FW_PH_dual_objective_output)
                            plot(collect(1:length(FW_PH_dual_objective_output[i])), FW_PH_dual_objective_output[i], xlims = (0,FW_PH_max_iter_used), xlabel = "iterations", ylabel = "Dual value resulting from FW-PH root node", label  = "FW-PH", legend = :bottomleft, title = "$(s)_$(i_fs_var)_$(i_ss_var)_$(i_const_num)_$(p)_")
                            if RNMDT_global_solved[i]
                                plot!(collect(1:FW_PH_max_iter_used), RNMDT_global_objective[i] .* ones(FW_PH_max_iter_used), label = "RNMDT")
                            end
                            savefig(output_link*"$(s)_$(i_fs_var)_$(i_ss_var)_$(i_const_num)_$(p)_" * " FW_PH_dual_value " * string(Dates.now()) * ".png")
                        end


                    end # for
                end # for
            end # for
        end # for

    # if we only solve the full scale problems using global solver
    elseif experiments_methods[1] == 1

        # define corresopnding data frame for the output dataset
        output_df = DataFrame( num_of_scen = Int[], num_fs_var = Int[], num_ss_var = Int[], num_const = Int[], p_RNMDT = Int[], primal_f = Float64[], primal_x = String[], primal_gap = Float64[])

        for s in n_scenarios
            for i_fs_var in fs_var
                for i_ss_var in ss_var
                    for i_const_num in const_num

                        # defining general intial parameters for the optimisation problem
                        initial_parameters = initialisation(s,i_fs_var,i_ss_var,i_const_num, p)

                        # generating the structure containing the constraints and objective related parameters
                        generated_parameters = parameters_generation(initial_parameters)

                        # solving the full scale problem using global solver
                        primal_problem = MIP_generation(initial_parameters, generated_parameters)
                        @suppress optimize!(primal_problem)
                        primal_problem_f = objective_value(primal_problem)
                        primal_problem_x = string(value.(primal_problem[:x][:,1]))
                        primal_problem_optimality_gap = MOI.get(primal_problem, MOI.RelativeGap())

                        push!(output_df, (s, i_fs_var, i_ss_var, i_const_num, 0 , primal_problem_f, primal_problem_x, primal_problem_optimality_gap))

                            # saving the updated database after solving each instance (in case there are some technical issues stopping the code)
                        XLSX.writetable(output_link*"experiments"*string(Dates.now())*".xlsx", output_df)

                    end # for
                end # for
            end # for
        end # for

    # if we only solve the RNMDT relaxations using global solver
    elseif experiments_methods[2] == 1

        # forming corresponding data frame
        output_df = DataFrame( num_of_scen = Int[], num_fs_var = Int[], num_ss_var = Int[], num_const = Int[], p_RNMDT = Int[], RNMDT_UB = Float64[], RNMDT_x = String[], RNMDT_time = Float64[], RNMDT_wy_gap = Float64[])

        for s in n_scenarios
            for i_fs_var in fs_var
                for i_ss_var in ss_var
                    for i_const_num in const_num

                        # defining starting value of the precision factor
                        p_fixed_value ? p = p_min_value : p = -1

                        # defining general intial parameters for the optimisation problem
                        initial_parameters = initialisation(s,i_fs_var,i_ss_var,i_const_num, p)

                        # generating the structure containing the constraints and objective related parameters
                        generated_parameters = parameters_generation(initial_parameters)

                        # defining the variables corresponing to the total time spend by global solver
                        # for solving the RNMDT isntances with different values of the precision factor p
                        rnmdt_g_time = 0

                        # while we did not exceed the total time allowed or did not reach the minimum value of the precision factor p
                        while (rnmdt_g_time < g_time_limit) && (p >= p_min_value)

                            RNMDT_relaxation = RNMDT_based_problem_generation(initial_parameters, generated_parameters)

                            # solving the RNMDT instance using the global solver
                            rndmt_init_time = time()
                            @suppress optimize!(RNMDT_relaxation)
                            rndmt_final_time = time() - rndmt_init_time
                            # updating the total time spent by the global solver to solve RNMDT instances
                            rnmdt_g_time += rndmt_final_time


                            # adding the line to the data base
                            push!(output_df, (s, i_fs_var, i_ss_var, i_const_num, p, objective_value(RNMDT_relaxation), string(value.(RNMDT_relaxation[:x][:,1])), rndmt_final_time, RNMDT_gap_computation( value.(RNMDT_relaxation[:y]), value.(RNMDT_relaxation[:w_RNMDT]))) )

                            # decreasing the value of the precision factor
                            p -= 1

                            # defining general intial parameters for the optimisation problem
                            initial_parameters = initialisation(s,i_fs_var,i_ss_var,i_const_num, p)

                            # generating the structure containing the constraints and objective related parameters
                            generated_parameters = parameters_generation(initial_parameters)

                            # saving the updated database after solving each instance (in case there are some technical issues stopping the code)
                            XLSX.writetable(output_link*"experiments"*string(Dates.now())*".xlsx", output_df)

                        end # while
                    end # for
                end # for
            end # for
        end # for

    # if we only solve the RNMDT relaxations using developed BnB_WF-PH method
    elseif experiments_methods[3] == 1

        # forming corresponding data frame
        output_df = DataFrame( num_of_scen = Int[], num_fs_var = Int[], num_ss_var = Int[], num_const = Int[], p_RNMDT = Int[], BnB_UB = Float64[], BnB_LB = Float64[], BnB_x = String[], BnB_time = Float64[], BnB_wy_gap = Float64[], BnB_nodes_explored = Int[] )


        for s in n_scenarios
            for i_fs_var in fs_var
                for i_ss_var in ss_var
                    for i_const_num in const_num

                        # defining starting value of the precision factor
                        p_fixed_value ? p = p_min_value : p = -1

                        # defining general intial parameters for the optimisation problem
                        initial_parameters = initialisation(s,i_fs_var,i_ss_var,i_const_num, p)

                        # generating the structure containing the constraints and objective related parameters
                        generated_parameters = parameters_generation(initial_parameters)

                        # defining the variables corresponing to the total time spent by the developed method
                        # for solving the RNMDT isntances with different values of the precision factor p
                        bnb_g_time = 0

                        # forming an array to collect the dual objective values resutling from iterations of FW-Ph
                        # applied to the instances with different values of the precision factor p
                        FW_PH_dual_objective_output = []

                        # the variable to keep the info on the max number of iterations used by FW-PH
                        # for the "fair" plot
                        FW_PH_max_iter_used = 0

                        # while we did not exceed the total time allowed or did not reach the minimum value of the precision factor p
                        while (bnb_g_time < g_time_limit) && (p >= p_min_value)


                            # if we need to plot the FW iterations for the root node
                            if initial_parameters.PH_SDM_parameters.PH_plot

                                # solving the RNMDT instance using the BnB+FW-PH method
                                bnb_p_init_time = time()
                                bnb_output, bnb_FW_iterations_output = bnb_solve(initial_parameters, non_ant_tol, tol_bb, integrality_tolerance)
                                bnb_p_final_time = time() - bnb_p_init_time

                                # collecting the data for plot
                                push!(FW_PH_dual_objective_output,bnb_FW_iterations_output[1])

                                # updating the max_number of the iteratiosn used by FW-PH for the different values of the precision factor p
                                length(bnb_FW_iterations_output[1])>FW_PH_max_iter_used ? FW_PH_max_iter_used = length(bnb_FW_iterations_output[1]) : true

                                # froming the dual_feasibility related output to be suitable for DataFrame
                                dual_feasibility_output = []
                                for i = 1:length(bnb_FW_iterations_output[2])
                                    push!(dual_feasibility_output, string(bnb_FW_iterations_output[2][i]))
                                end
                                # Printing the dual_deasibility_condtion and primal_dual_residual resulting from iterations of FW-PH applied to the root node
                                FW_PH_output = DataFrame( iteration = 1:length(bnb_FW_iterations_output[2]), dual_deasibility_cond = dual_feasibility_output, primal_dual_residual = bnb_FW_iterations_output[3])
                                XLSX.writetable(output_link*"$s : scen,$i_fs_var : fs_var, $i_ss_var : ss_var, $i_const_num : const, $p : prec_f" * string(Dates.now()) * "Dual feas cond and primal dual res " * ".xlsx", FW_PH_output)

                                # printing out the info on the identical etries in the feasibility set V
                                io1 = open(output_link*" $s : scen,$i_fs_var : fs_var, $i_ss_var : ss_var, $i_const_num : const, $p : prec_f" * string(Dates.now())* " feasibility set V identical entries count" * ".txt", "w")
                                for s = 1:length(bnb_FW_iterations_output[5])
                                        println(io1, "\n SCENARIO = $s : $(bnb_FW_iterations_output[5][s])")
                                end
                                close(io1)

                                # Printing the fesibility set V resulting from each iteration of FW-PH applied to the root node
                                io = open(output_link*" $s : scen,$i_fs_var : fs_var, $i_ss_var : ss_var, $i_const_num : const, $p : prec_f" * string(Dates.now())* " feasibility set V " * ".txt", "w")
                                    # printing only for the first scenario to ease representation
                                    for i = 1:length(bnb_FW_iterations_output[4])
                                        println(io, "\n ITERATION = $i : $(length(bnb_FW_iterations_output[4][i][1]))")
                                    end
                                    println(io, "--------------------------------------------------------")
                                    for i = 1:length(bnb_FW_iterations_output[4])
                                        println(io, "\n ITERATION = $i \n")
                                        println(io, bnb_FW_iterations_output[4][i][1])
                                        println(io,"\n")
                                    end
                                close(io)

                            else
                                # solving the RNMDT instance using the BnB+FW-PH method
                                bnb_p_init_time = time()
                                bnb_output = bnb_solve(initial_parameters, non_ant_tol, tol_bb, integrality_tolerance)
                                bnb_p_final_time = time() - bnb_p_init_time
                            end


                            bnb_g_time += bnb_p_final_time

                            # adding the line to the data base
                            push!(output_df, (s, i_fs_var, i_ss_var, i_const_num, p, bnb_output.UBDg, bnb_output.LBDg, string(bnb_output.soln_val), bnb_p_final_time, bnb_output.RNMDT_gap_wy, bnb_output.nodes_used))

                            # decreasing the value of the precision factor
                            p -= 1

                            # defining general intial parameters for the optimisation problem
                            initial_parameters = initialisation(s,i_fs_var,i_ss_var,i_const_num, p)

                            # generating the structure containing the constraints and objective related parameters
                            generated_parameters = parameters_generation(initial_parameters)

                            # saving the updated database after solving each instance (in case there are some technical issues stopping the code)
                            XLSX.writetable(output_link*"experiments"*string(Dates.now())*".xlsx", output_df)

                        end # while

                        # ploting the FW-PH dual objective related output
                        for i = 1:length(FW_PH_dual_objective_output)
                            plot(collect(1:length(FW_PH_dual_objective_output[i])), FW_PH_dual_objective_output[i], xlims = (0,FW_PH_max_iter_used), xlabel = "iterations", ylabel = "Dual value resulting from FW-PH root node", label  = "FW-PH", legend = :bottomleft, title = "$s : scen,$i_fs_var : fs_var, $i_ss_var : ss_var, $i_const_num : const, $(p_fixed_value ? p_min_value : -i) : prec_f")
                            savefig(output_link*"$s : scen,$i_fs_var : fs_var, $i_ss_var : ss_var, $i_const_num : const, $(p_fixed_value ? p_min_value : -i) : prec_f" * " FW_PH_dual_value " * string(Dates.now()) * ".png")
                        end

                    end # for
                end # for
            end # for
        end # for


    end # if

end

experiments_function(n_scenarios, fs_var, ss_var, const_num, experiments_methods, g_time_limit, p_value, output_link)

#link =  "/Users/nikitabelyak/Dropbox (Aalto)/branch-and-bound-caroe-and-schultz/"
#output_link =  "/scratch/work/belyakn1/BnB_p_lagrangian/"
#XLSX.writetable(output_link*"experiments"*string(Dates.now())*".xlsx", output_df)
