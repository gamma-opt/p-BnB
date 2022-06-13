"""
    RNMDT_gap_computation(y::Array{Float64}, w::Array{Float64})
Auxiliary function to compute the RNMDT-related difference
between y_i y_j and w_ij

"""
function RNMDT_gap_computation(y::Array{Float64}, w::Array{Float64})
    gap = 0
    for i = 1:size(y,1)
        for j = 1:size(y,1)
            for s = 1:size(y,2)
                gap += abs(y[i,s]*y[j,s] - w[i,j,s])

            end
        end
    end

    return gap

end
