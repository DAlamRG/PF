# This file contains the code neccesary to perform the pull-moves in a triangular lattice.



using LinearAlgebra


# First, I declare tha adjacency matrix for the protein configuration. 
# All I need is the protein length.
"""
    adjacencyM(n)

Given the protein´s length `n`; returns the adjacency matrix corresponding to the 
path formed by the protein.
"""
function adjacencyM(n)
    dv=zeros(Int8,n) # The diagonal contains zeros.
    ev=ones(Int8,n-1) # Uper diagonal contains only ones.
    M=Bidiagonal(dv, ev, :U)

    return M
end








