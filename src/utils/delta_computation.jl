"""
    delta_computation(X::Array{Float64})
returns parameter delta from the Caroe and Schultz algorithm
for the predefined values of the decision variables x

"""

function delta_computation(X::Array{Float64})
    delta = Array{Float64}(undef, size(X,1))
    [delta[i] = maximum(X[i,:]) - minimum(X[i,:]) for i = 1:size(X,1)]
end
