# This file contains the code neccesary to perform the pull-moves in a triangular lattice.


include("/Users/pedroruiz/Desktop/Diego/PF/HP2Dsquare/visualizeHP2D.jl")
include("/Users/pedroruiz/Desktop/Diego/PF/HP2Dsquare/HP-model2D-pullmoves.jl")
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





















# Now I write a function which checks for the validity of a configuration, which also depends on the geometry of the configuration.
"""
    validConfGeneral(N,ind,edo,HPlist,dir,geometry)

Given a 2D array of side `N`, a matrix encoding the aminoacids positions `edo`, an index `ind` , an array containing 
the sequence of H,P aminoacids `HPlist`, a direction `dir`, and a geometry; determines whether the protein structure is valid or not.
"""
function validConf2DGeneral(N,ind,edo,HPlist,dir,geometry)
    
    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)

    # Next, I iterate over the protein´s vertices, computing the "distance" to the next vertex. The direction in which I check the
    # structure is determined by `dir`.
    ans=true

    if geometry ==  square2D # Check for validty in square geometry.
        if dir == backwards
            if ind !=1
                for j in ind:-1:2
                    x1,y1=edo[j,:]
                    x1,y1=periodicInd2D(red,[x1,y1]) # Make the indices periodic.
                    x2,y2=edo[j-1,:]
                    x2,y2=periodicInd2D(red,[x2,y2]) # Make the indices periodic.
                    dx=abs(x1-x2)
                    dy=abs(y1-y2)
        
                    # Before assesing the validity of the sequence, we need to consider the case where the two elements being compared
                    # are located at the ends of the array `red`. In such case, the distance between them doesn´t necessarily exceed the
                    # limit.
        
                    xcond1= x1 == 1 && x2 == N
                    xcond2= x1 == N && x2 == 1
                    xcond= xcond1 || xcond2 
                    
                    ycond1= y1 == 1 && y2 == N
                    ycond2= y1 == N && y2 == 1
                    ycond= ycond1 || ycond2
        
    
                    if dx ==1 && dy ==1 
                        ans=false
                        break
                    
                    elseif dx > 1
                        if xcond == false
                            ans=false
                            break
                        end
                    end
                    if dy > 1
                        if ycond == false
                            ans=false
                            break
                        end
                    end
                end
            end
    
    
    
        elseif dir == forwards
            if ind != length(HPlist)
                for j in ind:length(HPlist)-1
                    x1,y1=edo[j,:]
                    x1,y1=periodicInd2D(red,[x1,y1]) # Make the indices periodic.
                    x2,y2=edo[j+1,:]
                    x2,y2=periodicInd2D(red,[x2,y2]) # Make the indices periodic.
                    dx=abs(x1-x2)
                    dy=abs(y1-y2)
        
                    # Before assesing the validity of the sequence, we need to consider the case where the two elements being compared
                    # are located at the ends of the array `red`. In such case, the distance between them doesn´t necessarily exceed the
                    # limit.
        
                    xcond1= x1 == 1 && x2 == N
                    xcond2= x1 == N && x2 == 1
                    xcond= xcond1 || xcond2 
                    
                    ycond1= y1 == 1 && y2 == N
                    ycond2= y1 == N && y2 == 1
                    ycond= ycond1 || ycond2
        
    
                    if dx ==1 && dy ==1 
                        ans=false
                        break
    
                    if dx > 1
                        if xcond == false
                            ans=false
                            break
                        end
                    end
                    if dy > 1
                        if ycond == false
                            ans=false
                            break
                        end
                    end
                end
            end
        end
    
    
    
    
    elseif geometry == triangular2D # Check the validty in triangular geometry.
        if dir == backwards
            if ind !=1
                for j in ind:-1:2
                    x1,y1=edo[j,:]
                    x1,y1=periodicInd2D(red,[x1,y1]) # Make the indices periodic.
                    x2,y2=edo[j-1,:]
                    x2,y2=periodicInd2D(red,[x2,y2]) # Make the indices periodic.
                    dx=abs(x1-x2)
                    dy=abs(y1-y2)
                    
                    # Before assesing the validity of the sequence, we need to consider the case where the two elements being compared
                    # are located at the ends of the array `red`. In such case, the distance between them doesn´t necessarily exceed the
                    # limit.
        
                    xcond1= x1 == 1 && x2 == N
                    xcond2= x1 == N && x2 == 1
                    xcond= xcond1 || xcond2 
                    
                    ycond1= y1 == 1 && y2 == N
                    ycond2= y1 == N && y2 == 1
                    ycond= ycond1 || ycond2
        
                    if dx > 1
                        if xcond == false
                            ans=false
                            break
                        end
                    end
                    if dy > 1
                        if ycond == false
                            ans=false
                            break
                        end
                    end
                end
            end
    
    
    
        elseif dir == forwards
            if ind != length(HPlist)
                for j in ind:length(HPlist)-1
                    x1,y1=edo[j,:]
                    x1,y1=periodicInd2D(red,[x1,y1]) # Make the indices periodic.
                    x2,y2=edo[j+1,:]
                    x2,y2=periodicInd2D(red,[x2,y2]) # Make the indices periodic.
                    dx=abs(x1-x2)
                    dy=abs(y1-y2)
        
                    # Before assesing the validity of the sequence, we need to consider the case where the two elements being compared
                    # are located at the ends of the array `red`. In such case, the distance between them doesn´t necessarily exceed the
                    # limit.
        
                    xcond1= x1 == 1 && x2 == N
                    xcond2= x1 == N && x2 == 1
                    xcond= xcond1 || xcond2 
                    
                    ycond1= y1 == 1 && y2 == N
                    ycond2= y1 == N && y2 == 1
                    ycond= ycond1 || ycond2
        
                    if dx > 1
                        if xcond == false
                            ans=false
                            break
                        end
                    end
                    if dy > 1
                        if ycond == false
                            ans=false
                            break
                        end
                    end
                end
            end
        end

    end

    return ans
end

























# Next, generalize the concept of nearest neighbors, according to the geometry.
"""
    nearestNeighbors2D(red,inds,geometry)

Given a 2D array `red`, and a couple of indices `inds`, and a geometry `geometry`; returns the values of the nearest neighbors to the given position,
in clockwise order. 
"""
function nearestNeighbors2D(red,inds,geometry)
    A=red
    x,y=inds
    x,y=periodicInd2D(red,[x,y])

    if geometry == square2D # Nearest neighbors in clockwise order.
        nn=Int8[periodicArr2D(A,[x-1,y]),periodicArr2D(A,[x,y+1]),periodicArr2D(A,[x+1,y]),periodicArr2D(A,[x,y-1])]
    elseif geomtry == triangular2D # Nearest neighbors in clockwise order.
        nn=Int8[periodicArr2D(A,[x-1,y]),periodicArr2D(A,[x,y+1]),periodicArr2D(A,[x+1,y+1]),periodicArr2D(A,[x+1,y])
        ,periodicArr2D(A,[x,y-1]),periodicArr2D(A,[x-1,y-1])]
    end
    
    return nn
end





















"""
    nearestNeighborsCoords2D(red,inds,geometry)

Given a 2D array `red`, and a couple of indices `inds`, and a geometry `geometry`; returns the coordinates of the nearest neighbors to 
the given position,in clockwise order.
"""
function nearestNeighborsCoords2D(red,inds,geometry)
    A=red
    x,y=inds
    x,y=periodicInd2D(red,[x,y])

    if geometry == square2D # Nearest neighbors in clockwise order.
        nnc=Int8[periodicInd2D(A,[x-1,y]),periodicInd2D(A,[x,y+1]),periodicInd2D(A,[x+1,y]),periodicInd2D(A,[x,y-1])]
    elseif geomtry == triangular2D # Nearest neighbors in clockwise order.
        nnc=Int8[periodicInd2D(A,[x-1,y]),periodicInd2D(A,[x,y+1]),periodicInd2D(A,[x+1,y+1]),periodicInd2D(A,[x+1,y])
        ,periodicInd2D(A,[x,y-1]),periodicInd2D(A,[x-1,y-1])]
    end
    
    return nnc
end























@enum NeighborsΔ::Int8 begin # Define nearest neighbors as enum type.
        upright=1
        right=2
        downright=3
        downleft=4
        left=5
        upleft=6
    end
"""
    countFirst2DΔ(N,edo,HPlist) 

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`,an array containing the sequence of H,P aminoacids `HPlist`, and
a geomtery; counts the number of possible pull moves for the first amino acid, 
recording the possible positions (for the triangle lattice).
"""
function countFirst2DGeneral(N,edo,HPlist)
    
    
    red=makeLattice(N,edo,HPlist)
    ind=1
    xplus,yplus=edo[ind+1,:]
    x,y=edo[ind,:]

    rel1= periodicInd2D(red,[x-1,y]) == periodicInd2D(red,[xplus,yplus])
    rel2= periodicInd2D(red,[x,y+1]) == periodicInd2D(red,[xplus,yplus])
    rel3= periodicInd2D(red,[x+1,y]) == periodicInd2D(red,[xplus,yplus])
    rel4= periodicInd2D(red,[x,y-1]) == periodicInd2D(red,[xplus,yplus]) # Conditions single out possible adjacent spaces to move the second amino acid to.

    vs=makecross2D(red,[x,y]) # These are the values of the nearest neighbors to `ind`.
    vscoords=[periodicInd2D(red,[x-1,y]),periodicInd2D(red,[x,y+1]),periodicInd2D(red,[x+1,y]),periodicInd2D(red,[x,y-1])] # coordinates
    # of nearest neighbors
    
    crossv1=makecross2D(red,[x-1,y])
    crossv2=makecross2D(red,[x,y+1])
    crossv3=makecross2D(red,[x+1,y])
    crossv4=makecross2D(red,[x,y-1]) # Nearest neighbors of the nearest neighbors.

    indices1=Int8[] # Here I will store the coordinates for `ind+1`.
    indices2=Neighbors[] # Here I store the possible coordinates for `ind`.
    if rel1 
        if vs[2] == 0 # `ind+1` might move into `vs[2]`.
            for k in [right,down]
                if crossv2[Int(k)] == 0
                    push!(indices1,2)
                    push!(indices2,k)
                end
            end
        end
        if vs[3] == 0
            for k in [right,down,left]
                if crossv3[Int(k)] == 0
                    push!(indices1,3)
                    push!(indices2,k)
                end
            end
        end
        if vs[4] == 0
            for k in [down,left]
                if crossv4[Int(k)] == 0
                    push!(indices1,4)
                    push!(indices2,k)
                end
            end
        end
    elseif rel2 
        if vs[1] == 0 
            for k in [up,left]
                if crossv1[Int(k)] == 0
                    push!(indices1,1)
                    push!(indices2,k)
                end
            end
        end
        if vs[3] == 0
            for k in [down,left]
                if crossv3[Int(k)] == 0
                    push!(indices1,3)
                    push!(indices2,k)
                end
            end
        end
        if vs[4] == 0
            for k in [up,down,left]
                if crossv4[Int(k)] == 0
                    push!(indices1,4)
                    push!(indices2,k)
                end
            end
        end
    elseif rel3
        if vs[1] == 0 
            for k in [up,right,left]
                if crossv1[Int(k)] == 0
                    push!(indices1,1)
                    push!(indices2,k)
                end
            end
        end
        if vs[2] == 0
            for k in [up,right]
                if crossv2[Int(k)] == 0
                    push!(indices1,2)
                    push!(indices2,k)
                end
            end
        end
        if vs[4] == 0
            for k in [up,left]
                if crossv4[Int(k)] == 0
                    push!(indices1,4)
                    push!(indices2,k)
                end
            end
        end
    elseif rel4 
        if vs[1] == 0 
            for k in [up,right]
                if crossv1[Int(k)] == 0
                    push!(indices1,1)
                    push!(indices2,k)
                end
            end
        end
        if vs[2] == 0
            for k in [up,right,down]
                if crossv2[Int(k)] == 0
                    push!(indices1,2)
                    push!(indices2,k)
                end
            end
        end
        if vs[3] == 0
            for k in [right,down]
                if crossv3[Int(k)] == 0
                    push!(indices1,3)
                    push!(indices2,k)
                end
            end
        end
    end
    
    # I now have two arrays with the information of where it is possible to move the first amino acids. First, I convert 
    # the numbers or enum types to coordinates. Then, making use of `validConfEnd2D` I make sure that the end amino acid 
    # `ind` doesn´t end up in the vicinity of another amino acid besides `ind+1`.

    coordinates1=[]
    coordinates2=[]

    for el in indices1
        push!(coordinates1,vscoords[el])
    end

    for k in 1:length(indices2)
        el=indices2[k]
        xv,yv=coordinates1[k]
        if el == up
            push!(coordinates2,periodicInd2D(red,[xv-1,yv]))
        elseif el == right
            push!(coordinates2,periodicInd2D(red,[xv,yv+1]))
        elseif el == down
            push!(coordinates2,periodicInd2D(red,[xv+1,yv]))
        elseif el == left
            push!(coordinates2,periodicInd2D(red,[xv,yv-1]))
        end
    end


    numpull=length(coordinates1)
    # Now we have the possible coordinates for the first two monomers, as well as the number of possible pull moves for the first 
    # amino acid. Fo easier use, I turn the arrays `coordinates1,coordinates2` into a `numpull×2` matrices.

    coordinates1=transpose(hcat(coordinates1...))
    coordinates2=transpose(hcat(coordinates2...))
    
    

    return (numpull,coordinates1,coordinates2)
end

















"""
    sharedNeighbors2D(red,inds,indsp)

Given a 2D array `red`, a couple of indices `inds` and `indsp`; returns the possible spaces to which `ind` might be moved to for 
the triangular geometry. 
"""
function sharedNeighbors2D(red,inds,indsp)
    sharedN=ones(Int8,6) # Possible shared neighbor spaces will have a value of zero.
    
    nnc=nearestNeighborsCoords2D(red,inds,triangular2D) # Nearest neigbors to our index `ind`.
    nn=nearestNeighbors2D(red,inds,triangular2D)
    
    rel1= nnc[1] == periodicInd2D(red,indsp)
    rel2= nnc[2] == periodicInd2D(red,indsp)
    rel3= nnc[3] == periodicInd2D(red,indsp)
    rel4= nnc[4] == periodicInd2D(red,indsp)
    rel5= nnc[5] == periodicInd2D(red,indsp)
    rel6= nnc[6] == periodicInd2D(red,indsp) # Possible relative positions between `inds,indsp`.

    if rel1 # Fill `diagspaces` according to the relative positions.
        sharedN[2]= nn[2]
        sharedN[6]= nn[6]
    elseif rel2 
        sharedN[1]= nn[1]
        sharedN[3]= nn[3]
    elseif rel3 
        sharedN[2]= nn[2]
        sharedN[4]= nn[4]
    elseif rel4 
        sharedN[3]= nn[3]
        sharedN[5]= nn[5]
    elseif rel5
        sharedN[4]= nn[4]
        sharedN[6]= nn[6]
    elseif rel6
        sharedN[5]= nn[5]
        sharedN[1]= nn[1]
    end

    return sharedN
end















"""
    countMiddle2DΔ(N,ind,edo,HPlist) 

Given a 2D array size `N`, an index `ind`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; counts the number of possible pull moves for the middle amino acids, and 
stores the possible coordinates (for triangle lattice).
"""
function countMiddle2DΔ(N,ind,edo,HPlist)
    
    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)
    

    # Now we check wether a pull-move is possible for the selected vertex. A shared neighbor space between `ind` and either `ind+1` or `ind-1`
    # must be empty.
    # First I check where `ind+1,ind-1` are.
    xplus,yplus=edo[ind+1,:] # Position of `ind+1`
    xminus,yminus=edo[ind-1,:] # Position of `ind-1`
    x,y=edo[ind,:] # Position of `ind`


    # Next we store the possible free spaces for `ind` when pulling backwards.
    freespaces1=sharedNeighbors2D(red,[x,y],[xplus,yplus])
    
    
    # Store possible free spaces when pulling forwards.
    freespaces2=sharedNeighbors2D(red,[x,y],[xminus,yminus])

    nnc=nearestNeighborsCoords2D(red,[x,y],triangular2D)
    

    # Now I write the possible coordinates for `ind`.
    coordinates1=[]
    for k in [upright,right,downright,downleft,left,upleft]
        space=freespaces1[Int(k)]
        if space == 0 
            push!(coordinates1,nnc[Int(k)])
        else
            continue
        end
    end

    coordinatesi=[]
    for k in [upright,right,downright,downleft,left,upleft]
        space=freespaces2[Int(k)]
        if space == 0 
            push!(coordinatesi,nnc[Int(k)])
        else
            continue
        end
    end
    


    # Now I have in `coordinates1,coordinatesi` the potential new positions for `ind`. The length of these arrays is 
    # the number of possible pull moves por the given index.
    numpull1=length(coordinates1)
    numpull2=length(coordinatesi)
    numpull=numpull1+numpull2

    # I turn the 1 dimensional arrays into matrices for easier use.
    coordinates1=transpose(hcat(coordinates1...))
    coordinatesi=transpose(hcat(coordinatesi...))

    return (numpull,coordinates1,coordinatesi) # Returns everything neccesary to reproduce the final configuration.
end









