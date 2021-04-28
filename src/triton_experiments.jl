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

## Constructing the experiments
# the structure that will collect the experiments results
output_df = DataFrame( num_of_scen = Int[], num_fs_var = Int[], num_ss_var = Int[], num_const = Int[], p_RNMDT = Int[], primal_f = Float64[], primal_x = String[], primal_gap = Float64[], RNMDT_UB = Float64[], RNMDT_x = String[], RNMDT_time = Float64[], RNMDT_wy_gap = Float64[], BnB_UB = Float64[], BnB_LB = Float64[], BnB_x = String[], BnB_time = Float64[], BnB_wy_gap = Float64[], BnB_nodes_explored = Int[] )

scenarios = [5,10,15]
#scenarios = [15]
fs_var = [5, 7, 10]
#fs_var = [5]
output_link =  "/scratch/work/belyakn1/BnB_p_lagrangian/"

for s in scenarios
    for i_fs_var in fs_var


        p = -1
        # defining general intial parameters for the optimisation problem
        initial_parameters = initialisation(s,i_fs_var,i_fs_var,i_fs_var, p)

        # generating the structure containing the constraints and objective related parameters
        generated_parameters = parameters_generation(initial_parameters)

        primal_problem = MIP_generation(initial_parameters, generated_parameters)
        optimize!(primal_problem)
        primal_problem_f = objective_value(primal_problem)
        primal_problem_x = string(value.(primal_problem[:x][:,1]))
        primal_problem_optimality_gap = MOI.get(primal_problem, MOI.RelativeGap())

        bnb_g_time = 0
        rnmdt_g_time = 0

            while (bnb_g_time < 3600) && (p >= -12)

                RNMDT_relaxation = RNMDT_based_problem_generation(initial_parameters, generated_parameters)

                if rnmdt_g_time < 3600
                    rndmt_init_time = time()
                    optimize!(RNMDT_relaxation)
                    rndmt_final_time = time() - rndmt_init_time
                    rnmdt_g_time += rndmt_final_time
                    rndmt_happened = true
                else
                    rndmt_happened = false
                end






                bnb_p_init_time = time()
                bnb_output = bnb_solve(initial_parameters, non_ant_tol, tol_bb, integrality_tolerance)
                bnb_p_final_time = time() - bnb_p_init_time

                bnb_g_time += bnb_p_final_time

                if rndmt_happened
                    push!(output_df, (s, i_fs_var, i_fs_var, i_fs_var, p, primal_problem_f, primal_problem_x, primal_problem_optimality_gap, objective_value(RNMDT_relaxation), string(value.(RNMDT_relaxation[:x][:,1])), rndmt_final_time, RNMDT_gap_computation( value.(RNMDT_relaxation[:y]), value.(RNMDT_relaxation[:w_RNMDT])), bnb_output.UBDg, bnb_output.LBDg, string(bnb_output.soln_val), bnb_p_final_time, bnb_output.RNMDT_gap_wy, bnb_output.nodes_used))
                else
                    push!(output_df, (s, i_fs_var, i_fs_var, i_fs_var, p, primal_problem_f, primal_problem_x, primal_problem_optimality_gap, 0.0, "NaN", 0.0, 0.0, bnb_output.UBDg, bnb_output.LBDg, string(bnb_output.soln_val), bnb_p_final_time,  bnb_output.RNMDT_gap_wy, bnb_output.nodes_used ))
                end
                p -= 1

                # defining general intial parameters for the optimisation problem
                initial_parameters = initialisation(s,i_fs_var,i_fs_var,i_fs_var,p)

                # generating the structure containing the constraints and objective related parameters
                generated_parameters = parameters_generation(initial_parameters)

                XLSX.writetable(output_link*"experiments"*string(Dates.now())*".xlsx", output_df)
            end

    end

end
#link =  "/Users/nikitabelyak/Dropbox (Aalto)/branch-and-bound-caroe-and-schultz/"
output_link =  "/scratch/work/belyakn1/BnB_p_lagrangian/"
XLSX.writetable(output_link*"experiments"*string(Dates.now())*".xlsx", output_df)
