## This file contains the fucntion required to extract the infomation about 
## the constraints and objective functions coefficents from the pooling problem

function polling_problem_parameters_extraction(model::JuMP.Model, bub)
    #const MOI = MathOptInterface

    # collecting the list of the variables 
    vars = JuMP.all_variables(model)
    # collecting the list of indicies
    indices = MOI.get(model, MOI.ListOfVariableIndices())
    # mapping the Variablesref to the index
    Variables = Array{pooling_variable}(undef, length(vars))
    for i = 1:length(vars)
        Variables[i] = pooling_variable(string(vars[i]), indices[i].value)
    end



    # extracrting the data for linear objective
    # this only works if we have a linear objective function (the model has a ScalarAffineFunction as its objective)
    obj = MOI.get(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    Objective_c = zeros(length(Variables))
    for term in obj.terms
        Objective_c[term.variable.value] = term.coefficient
    end

    # extracting the data from constraints
    # list all the constraints present in the model
    cons = MOI.get(model, MOI.ListOfConstraintTypesPresent())
    #@show(cons)

    Quadratic_constraints = []
    Box_constraints = []
    Affine_constraints = Array{af_const}(undef, 2)
 

    # go through all the types of constraints 
    for i = 1:length(cons)
        
        # get the constraint indices for this combination of F(unction) in S(et)
        F = cons[i][1]
        S = cons[i][2]

        # We get two constraint indices (stored in the array ci),
        ci = MOI.get(model, MOI.ListOfConstraintIndices{F,S}())   


        # if the constraint is quadratic 
        if F == MathOptInterface.ScalarQuadraticFunction{Float64}
        # create axuliairy matrix inside the loop as only here we know how many quadratic constraints we have 
        Q_c = Array{quad_const}(undef, length(ci))

            # got through all the constraints of this type
            for j = 1:length(ci)
                cij = ci[j]

                # to get the function and set corresponding to this constraint (index):
                moi_backend = JuMP.backend(model)
                f = MOI.get(moi_backend, MOI.ConstraintFunction(), cij)
                
                # create a matrix to store the quadratic coefficeints of constraint j
                Qj = zeros(length(vars),length(vars))
                # create a vector to store the affine coefficeints of constraint j
                Aj = zeros(length(vars))

                # assemble the quadratic coefficeints matrix
                for term in f.quadratic_terms
                    Qj[term.variable_1.value, term.variable_2.value] = term.coefficient
                end

                # assemble the affine coefficeints vector
                for term in f.affine_terms
                    Aj[term.variable_index.value] = term.coefficient
                end

                # remember the constant 
                Cj = f.constant

                #create correspondent strucutre 
                quadratic_cosntraint = quad_const(Qj, Aj, string(S), Cj)
                # add the structure
                Q_c[j] =  quadratic_cosntraint
            end
        # Save the inner loop variable Q_c values to in the global variable Quadratic_constraints
        Quadratic_constraints = Q_c
        elseif F == MathOptInterface.ScalarAffineFunction{Float64} # if the constraint is afine
            # create matrix to keep all linear coeffiecents
            A = Array{Float64}(undef, length(ci), length(vars))
            # create vector of the RHS values 
            B = Array{Float64}(undef, length(ci))
            for j = 1:length(ci)
                cij = ci[j]
                # to get the function and set corresponding to this constraint (index):
                moi_backend = JuMP.backend(model)
                f = MOI.get(moi_backend, MOI.ConstraintFunction(), cij)
                # create auxiliary vector that would be j-row of matrix A
                aj = zeros(length(vars))
                for term in f.terms
                    aj[term.variable.value] = term.coefficient
                end
                
                # To get the corresponding first entry bj of b = [b1; ...; bm] vector, 
                # we have to look at the constraint set of that same constraint index cij:
                s = MOI.get(moi_backend, MOI.ConstraintSet(), cij)

                # store correspondet j-row and j value fo the coeeficent matricx A and RHS values vector B
                A[j,:] = aj
                B[j] = (S == MathOptInterface.EqualTo{Float64} ? s.value : s.upper)
                
                if S == MathOptInterface.EqualTo{Float64}
                    Affine_constraints[1] = af_const(A, string(S), B)
                else
                    Affine_constraints[2] = af_const(A, string(S), B)
                end

            end
            
        else # if it is half-box constraint
            # create axuliairy matrix inside the loop as only here we know how many quadratic constraints we have 
            B_c = Array{box_const}(undef, length(ci))

            for j = 1:length(ci)
                cij = ci[j]
                
                # to get the function and set corresponding to this constraint (index):
                moi_backend = JuMP.backend(model)
                f = MOI.get(moi_backend, MOI.ConstraintFunction(), cij)
                var_index = f.value

                # To get the corresponding first entry bj of b = [b1; ...; bm] vector, 
                # we have to look at the constraint set of that same constraint index cij:
                s = MOI.get(moi_backend, MOI.ConstraintSet(), cij)
                RHS = (S == MathOptInterface.LessThan{Float64} ? s.upper : s.lower)
                B_c[j] = box_const( var_index, string(S), RHS)
            end
            Box_constraints = B_c
        end
        
    end
   

    return Variables, Objective_c, Affine_constraints, Quadratic_constraints, Box_constraints, bub
end