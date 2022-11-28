# Write here the local repository address
#src_link  =  "/scratch/work/belyakn1/BnB_p_lagrangian/src/"
src_link  =  "/Users/nikitabelyak/Dropbox (Aalto)/BnB_Caroe-and-Schultz/src/"
#src_link =  "/Users/Fabricio/Documents/GitHub/BnB_Caroe-and-Schultz/src/"

#cd(chop(src_link, tail = 4))
using Pkg
Pkg.activate(".")
#Pkg.update()

Pkg.instantiate()

using JuMP, Gurobi, Random, LinearAlgebra, SparseArrays, Suppressor, Ipopt, Profile, Base.Threads
using DataFrames, XLSX, Dates, CSV, MathOptInterface

# set unique envinronment for Gurobi
const GRB_ENV = Gurobi.Env()

# check whether the folder for the "today" experiments exists
# and if not create one 
if !isdir(chop(src_link, tail = 4) * "experiments_" * string(Dates.today()))
    mkdir(chop(src_link, tail = 4) * "experiments_" * string(Dates.today()))
end

output_link = chop(src_link, tail = 4) * "experiments_" * string(Dates.today()) * "/"

## Defining the initial parameters for the experiments

# random seed 
r_seed = 0
# the number of scenarios
n_scenarios = [10]
# the number of the first stage variables / per scenario
fs_var = [5]
# the number of the second stage variabes / per scenario
ss_var = [20]
# the number of the constraints / per scenario
const_num = [5]
# Methods to be used [full scale, RNMDT, BnB]: 1 = use, 0 = don't use
experiments_methods = [1, 0, 0]
# The maximum time limit for the set of instances
g_time_limit = 3600
# precision factor values to be considered for RNMDT
p_value = [-1]
# Indicator whether bundle method is used (false) or Frank-Wolfe Progressive Hedging (true) 
AL = true
# Minmisation (true) or Maximisation(false) primal problem 
pp_min_max = false

include(src_link*"initialization.jl")


# Calling the experiments function
# (the results will be safed in the folder with correspoding date as a title)
experiments_function(n_scenarios, fs_var, ss_var, const_num, experiments_methods, g_time_limit, p_value, output_link, AL, r_seed)

