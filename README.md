
This repository contains the implememntation of the [Carøe and Schultz](https://www.sciencedirect.com/science/article/pii/S0167637798000509) propsed branch and bound method encased around [p-Lagrangian  relaxation](https://link.springer.com/article/10.1007/s10898-022-01138-y) method to prevent non-convergence due to mixed-integer nature of the input instances. 

# Usage 
1. ## Initialisation of the parameteres for the algorithms and parameters for the Random fucntion

2. ## Initialisation of the parameteres for the test instances 
    - ### 2.a. Randomly generated non-convex MIQCQP instances
        Inside the `src` directory we have `experiments.jl` file. In this file please define the parameters for the instance(s) that you want to test 
        ```
        ## Defining the initial parameters for the experiments

        # the number of scenarios
        n_scenarios = [20]
        # the number of the first stage variables / per scenario
        fs_var = [5]
        # the number of the second stage variabes / per scenario
        ss_var = [10]
        # the number of the constraints / per scenario
        const_num = [10]
        # Methods to be used [full scale, RNMDT, BnB]: 1 = use, 0 = don't use
        experiments_methods = [0, 0, 1]
        # The maximum time limit for the set of instances
        g_time_limit = 3600
        # precision factor values to be considered for RNMDT
        p_value = [-1]
        ```
    - ### 2.b. Pooling problem with one pool based stochastic instances
3. ## Running example instance
    Inside the `src` directory we have `experiments.jl` file that contains the example on how to run file using 