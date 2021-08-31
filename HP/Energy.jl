# This script computes the energy asociated with a protein structure.


include("./HP-pullmoves.jl")











"""
    interaction(amin1::Amin,amin2::Amin,pfmodel::PF_model)::Float64
Given a two aminoacids `amin1,amin2` and a model for the folding process `pfmodel` (models available 
are `HP1_model,HP2_model,HP3_model,HPNX_model,hHPNX_model,YhHX_model`); returns a value for the interaction between the given 
aminoacids.
"""
function interaction(amin1::Amin,amin2::Amin,pfmodel::PF_model)::Float64
    i = pfmodel.dict[amin1]
    j = pfmodel.dict[amin2]
    matrix = pfmodel.int_matrix
    return matrix[i,j]
end














"""
    energyNeighbors(red,edo,ind,HPlist,geometry)

Given a 2D/3D array `red`, the list of aminoacid poitions `edo` , an index on the chain `ind`, the list of aminoacid types `HPlist` and a 
geometry; returns the values for the topological neighbors of `ind` for the given geometry. These list of values excludes the values of 
the indices immediately adjacent to our index.
"""
function energyNeighbors(red,edo,ind,HPlist,geometry)
    nnc = nearestNeighborsCoords(red,edo[ind,:],geometry) # Coordiantes of nearest neigbors to our index `ind`.
    nn = nearestNeighbors(red,edo[ind,:],geometry) # Value of nearest neighbors to our index.

    if ind == 1
        aminvals = Amin[] # Amin values of topological neighbors of inds, which are not adjacent to it.
        for i in 1:length(nnc)
            if nnc[i] ≠ periodicInd(red,edo[ind+1,:],length(edo[1,:])) && nn[i] ≠ 0
                push!(aminvals,amin_dict[nn[i]]) 
            end
        end

        return aminvals

    elseif ind == length(HPlist)
        aminvals = Amin[] # Amin values of topological neighbors of inds, which are not adjacent to it.
        for i in 1:length(nnc)
            if nnc[i] ≠ periodicInd(red,edo[ind-1,:],length(edo[1,:])) && nn[i] ≠ 0
                push!(aminvals,amin_dict[nn[i]]) 
            end
        end

        return aminvals

    else
        aminvals = Amin[] # Amin values of topological neighbors of inds, which are not adjacent to it.
        for i in 1:length(nnc)
            if (nnc[i] ≠ periodicInd(red,edo[ind-1,:],length(edo[1,:]))) && (nnc[i] ≠ periodicInd(red,edo[ind+1,:],length(edo[1,:]))) && nn[i] ≠ 0
                push!(aminvals,amin_dict[nn[i]]) 
            end
        end

        return aminvals
    end
end














"""
    energy(N,edo,HPlist,geometry,pfmodel)

Given a 2D/3D array size `N`, a matrix encoding the aminoacids positions `edo`, an array containing the sequence of aminoacids `HPlist`,the 
geometry of the lattice and an interaction model `pfmodel`; returns the energy of the configuration.
"""
function energy(N,edo,HPlist,geometry,pfmodel::PF_model)
    red = makeLattice(N,edo,HPlist)
    
    # Next, I iterate over the protein´s vertices, computing the energy contribution by each of the non-adjacent bonds, according to the model.
    ε = 0

    for ind in 1:length(HPlist)
        pos = periodicArr(red,edo[ind,:],length(edo[1,:]))
        amin1 = amin_dict[pos]
        listamins = energyNeighbors(red,edo,ind,HPlist,geometry)
        for amin2 in listamins
            ε = ε + interaction(amin1,amin2,pfmodel)
        end
    end
    return ε/2
end





# Test the energy function
# edo_test = [[5 6];[5 7];[6 7];[6 8];[5 8];[5 9];[6 9];[7 9];[7 8];[7 7];[7 6];[6 6]]

# HPlist_test= [H,P,P,H,H,P,H,P,P,H,P,H]

# @show(energy(11,edo_test,HPlist_test,square2D,HP3_model))