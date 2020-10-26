"""
    gurobi_solver_parameters
Stores the values for the parameters of the Gurobi Optimizer. Has the following fields:

* `Method::Int`:                     Relates to the algorithm used to solve continuous models
                                     https://www.gurobi.com/documentation/9.0/refman/method.html#parameter:Method

* `OutputFlag::Int`:                 Relates to the solver output control
                                     https://www.gurobi.com/documentation/9.0/refman/outputflag.html#parameter:OutputFlag

* `MIPGap::Int`:                     Relates to the relative MIP optimality gap
                                     https://www.gurobi.com/documentation/9.0/refman/mipgap2.html#parameter:MIPGap

* `Threads::Int`:                    Relates to the number of parallel threads to use
                                     https://www.gurobi.com/documentation/9.0/refman/threads.html#parameter:Threads

* `NonConvex::Int`:                  Relates to the control how to deal with non-convex quadratic programs
                                     https://www.gurobi.com/documentation/9.0/refman/nonconvex.html#parameter:NonConvex

"""

mutable struct gurobi_solver_parameters
    Method::Int
    OutputFlag::Int
    MIPGap::Int
    Threads::Int
    NonConvex::Int
end
