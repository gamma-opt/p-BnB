
This repository contains the implememntation of the [Carøe and Schultz](https://www.sciencedirect.com/science/article/pii/S0167637798000509) propsed branch and bound method encased around [p-Lagrangian  relaxation](https://link.springer.com/article/10.1007/s10898-022-01138-y) method to prevent non-convergence due to mixed-integer nature of the input instances. 

# Usage 
1. ## Initialisation of the parameteres for the algorithms and parameters for the Random fucntion
    Inside the `src>utils` directory we have `initialisation_function.jl` file. Inside this file pleas fill in all the necessary parameters except for those denoted with `pooling section parameters`
    ```julia
    #-----------------------------------------------------------------------------
    # pooling problem parameters
    g_pool_problem_is_used = true

    pp_primal_pool_problem_link = chop(src_link, tail = 4) * "pooling_problem_data"
    pp_start_value_pooling_nodes = [1.0]
    g_flow_cost = 10
    g_pool_cost = 10

    # stochastic pooling problem parameters
    stoch_pp_pool_max_capacity = 10.0
    stoch_pp_flow_max_capacity = 10.0
    stoch_pp_demand_min_value = 0.0 
    stoch_pp_demand_max_value = 5.0
    stoch_pp_source_node_min_value = 1.0
    stoch_pp_source_node_max_value = 3.0

    g_stoc_ppm = stochastic_pooling_problem_parameters(stoch_pp_pool_max_capacity, stoch_pp_flow_max_capacity, stoch_pp_demand_min_value, stoch_pp_demand_max_value, stoch_pp_source_node_min_value, stoch_pp_source_node_max_value)

    g_pooling_problem_parameters = pooling_problem_parameters(pp_primal_pool_problem_link, pp_start_value_pooling_nodes, g_flow_cost, g_pool_cost, g_stoc_ppm)
    #-----------------------------------------------------------------------------
    ```
2. ## Initialisation of the parameteres for the test instances 
    - ### 2.a. Randomly generated non-convex MIQCQP instances
        - Inside the `src` directory we have `experiments.jl` file. In this file please define the parameters for the instance(s) that you want to test 
            ```julia
            ## Defining the initial parameters for the experiments

            # the number of scenarios
            n_scenarios = [20]
            # the number of the first stage variables / per scenario
            fs_var = [5]
            # the number of the second stage variabes / per scenario
            ss_var = [10]
            # the number of the constraints / per scenario
            const_num = [10]
            # Methods to be used [full scale, RNMDT, BBPH]: 1 = use, 0 = don't use
            experiments_methods = [0, 0, 1]
            # The maximum time limit for the set of instances
            g_time_limit = 3600
            # precision factor values to be considered for RNMDT
            p_value = [-1]
            ```
        - Inside the `src>utils` directory we have `initialisation_function.jl` file. Please set up the parameter `g_pool_problem_is_used = true`
    - ### 2.b. Pooling problem with one pools layer based stochastic instances
        If you want the test instance to be generated based on exisiting pooling problem with one pools layer then 
        - in the folder `pooling_problem_data` fill in the data about the pooling problem you are going to use regarding arcs and nodes in `arcs.csv` and `nodes.csv` files respectively.
        - go back to step 1. and fill in the data for the section denoted with `pooling section parameters`
        - do the step 2.a but consider that the parameters `fs_var`, `ss_var` and `const_num` will be altered to match the values from  resulting stochastic instance based on the pooling problem.
3. ## Running example instance
    Inside the `src` directory we have `experiments.jl` file that contains the example on how to run the algorithms using function 
    `experiments_function`. If you completeted step 2.a you can run it as 
    ```julia 
    experiments_function(n_scenarios, fs_var, ss_var, const_num, experiments_methods, g_time_limit, p_value, output_link)
    ```
    Otherwise, you will have to manually substitue all the values as for instance the following line 
    ```julia 
    experiments_function(2, 10, 10, 10, [0,0,1], g_time_limit, [-1], output_link)
    ```