"""
    gurobi_solver_parameters
Stores the values for the parameters of the Gurobi Optimizer. Has the following fields:

* `Method::Int`:                     Relates to the algorithm used to solve continuous models
                                     https://www.gurobi.com/documentation/9.0/refman/method.html#parameter:Method

* `OutputFlag::Int`:                 Relates to the solver output control
                                     https://www.gurobi.com/documentation/9.0/refman/outputflag.html#parameter:OutputFlag

* `IntFeasTol::Float64`:             Relates to the solver integrality tolerance
                                     https://www.gurobi.com/documentation/9.1/refman/tolerances_and_user_scalin.html

* `FeasibilityTol::Float64`:         Relates to the solver primal feasibility tolerance
                                     https://www.gurobi.com/documentation/9.1/refman/tolerances_and_user_scalin.html

* `OptimalityTol::Float64`:          Relates to the solver dual feasibility tolerance
                                     https://www.gurobi.com/documentation/9.1/refman/tolerances_and_user_scalin.html

* `Threads::Int`:                    Relates to the number of parallel threads to use
                                     https://www.gurobi.com/documentation/9.0/refman/threads.html#parameter:Threads

* `NonConvex::Int`:                  Relates to the control how to deal with non-convex quadratic programs
                                     https://www.gurobi.com/documentation/9.0/refman/nonconvex.html#parameter:NonConvex

* `NumericFocus::Int`:               Parameter controls how the solver manages numerical issues
                                     https://www.gurobi.com/documentation/9.1/refman/making_the_algorithm_less_.html

"""

mutable struct gurobi_solver_parameters
    Method::Int
    OutputFlag::Int
    IntFeasTol::Float64
    FeasibilityTol::Float64
    OptimalityTol::Float64
    Threads::Int
    NonConvex::Int
    NumericFocus::Int
end
