"""
    fucntion arrays_equal(x, y)
    The function that compares wether nxn dimenational arrays
    of Float64 elements x and y are equal and
    returns true if it is the case and false otherwise.
"""

function arrays_equal(x, y)

    # introducing an auxiliary array of boolean values
    # at which (i,j) elemnt corresponds to whether the following statement is true or false
    # x[i,j] = y[i,j]
    check_array = x .== y

    # reshaping the array into a vector
    check_array = reshape(check_array, size(check_array,1)*size(check_array,2))

    # checing whether all the elemnts of check_array) have "true" values
    return sum(check_array)/length(check_array) == 1

end


"""
    fucntion array_elements_equal( x, y )
    The function that compares wether all the elemnts of the
    1-dimentional array x and y are indentical and returns true
    if they are and false otherwise.
"""

function array_elements_equal( x, y )

    # defining the variable corresposndent to the final verdict
    elements_equal = true

    # claucliating the number of elements in x (same as in y)
    l = length(x)

    # going through all the elemnts
    for i = 1:length(x)

        # checking wether correspondent elemtns of the x and y coincide
        # and updating the value of elements_equal accordingly
        elements_equal = elements_equal & arrays_equal(x[i], y[i])

        # stopping the algorithm once we found problematic element
        if !elements_equal
            return elements_equal
        end

    end

    return elements_equal
end

"""
    function feasibility_set_contains_element(V, new_element)
    The function that checks the appearnce of the new_element
    in the vector V and returns true if it is there and false otherwise
"""

function feasibility_set_contains_element(V, new_element)

    # defining the variable corresponding to the final verdict
    already_exists = false

    # calculating the number of the elements in V
    i = length(V)

    # going through all the elemtns till we rich the end of them
    # or find the new_element amonth them
    while !already_exists &&  (i > 0)

        # checking whether the new_element is equal to the element V[i]
        already_exists = array_elements_equal(V[i], new_element)

        # stopping the algorithm if we face coincidence
        if already_exists
            return already_exists
        # moving to the next elemnt otherwise
        else
            i -= 1
        end
    end

    return already_exists
end


# Some testing

#t1 = [1.0  1.0; 1.0 1.0]
#t2 = [1.0 1.0; 2.0 1.0]

# arrays_equal(t1, t2)

# test1 = [t1, t2, [1.0 1.0]]
# test2 = [t1, t2, [1.0 2.0]]

# array_elements_equal(test1, test2)

# Vtest = [test1, test2]

# feasibility_set_contains_element(Vtest, test1)
