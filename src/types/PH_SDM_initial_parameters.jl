
"""
    PH_SDM_input
Stores the input parameters for Frank-Wolfe Progressive Hedging (FW_PH) and Simplical Decomposition Method (SDM). Has the following fields:

* `α::Float64`:                                                     Parameter that affects initial linearisation point of the SDM method
* `PH_tol::Float64`:                                                Tolerance used for stopping criterion in FW-PH
* `PH_max_iter`:                                                    Maximum number of iterations used in FW-PH
* `PH_plot`:                                                        The indicator whether we should or not plot the FW-PH porgrees for the root node of the BnB
* `PH_max_time`:                                                    The maximum time allowed per instance
* `SDM_tolerance`:                                                  Tolerance used for stopping criterion in SDM
* `SDM_max_iter`:                                                   Maximum number of iterations used in SDM

"""
mutable struct PH_SDM_input
    α::Float64
    PH_tol::Float64
    PH_max_iter::Int
    PH_plot::Bool
    PH_max_time::Float64
    SDM_tolerance::Float64
    SDM_max_iter::Int
end
