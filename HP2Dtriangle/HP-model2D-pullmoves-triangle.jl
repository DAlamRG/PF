# This file contains the code neccesary to perform the pull-moves in a triangular lattice.


include("/Users/pedroruiz/Desktop/Diego/PF/HP2Dsquare/visualizeHP2D.jl")
include("/Users/pedroruiz/Desktop/Diego/PF/HP2Dsquare/HP-model2D-pullmoves.jl")
using LinearAlgebra




















# First, I declare tha adjacency matrix for the protein configuration. 
# All I need is the protein length. (This function might still be needed).
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
    elseif geometry == triangular2D # Nearest neighbors in clockwise order.
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
    elseif geometry == triangular2D # Nearest neighbors in clockwise order.
        nnc=[periodicInd2D(A,[x-1,y]),periodicInd2D(A,[x,y+1]),periodicInd2D(A,[x+1,y+1]),periodicInd2D(A,[x+1,y])
        ,periodicInd2D(A,[x,y-1]),periodicInd2D(A,[x-1,y-1])]
    end
    
    return nnc
end



























"""
    sharedNeighbors2D(red,inds,indsp)

Given a 2D array `red`, a couple of indices `inds` and `indsp`; returns the empty shared spaces which `inds,indsp` share for 
the triangular geometry. 
"""
function sharedNeighbors2D(red,inds,indsp)
    sharedN=ones(Int8,6) # Empty shared neighbor spaces will have a value of zero.
    
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
    excludedNeighbors2D(red,inds,indsp)

Given a 2D array `red`, a couple of indices `inds` and `indsp`; returns the value of the nearest neigbors to `inds` that are not shared by `indsp`.
"""
function excludedNeighbors2D(red,inds,indsp)
    notsharedN=ones(Int8,6) # Empty not-neighbor spaces will have a value of zero.
    
    nnc=nearestNeighborsCoords2D(red,inds,triangular2D) # Nearest neigbors to our index `ind`.
    nn=nearestNeighbors2D(red,inds,triangular2D) # Value of closest neighbors.
    
    rel1= nnc[1] == periodicInd2D(red,indsp)
    rel2= nnc[2] == periodicInd2D(red,indsp)
    rel3= nnc[3] == periodicInd2D(red,indsp)
    rel4= nnc[4] == periodicInd2D(red,indsp)
    rel5= nnc[5] == periodicInd2D(red,indsp)
    rel6= nnc[6] == periodicInd2D(red,indsp) # Possible relative positions between `inds,indsp`.

    if rel1 # Fill `diagspaces` according to the relative positions.
        notsharedN[3]= nn[3]
        notsharedN[4]= nn[4]
        notsharedN[5]= nn[5]
    elseif rel2 
        notsharedN[4]= nn[4]
        notsharedN[5]= nn[5]
        notsharedN[6]= nn[6]
    elseif rel3 
        notsharedN[5]= nn[5]
        notsharedN[6]= nn[6]
        notsharedN[1]= nn[1]
    elseif rel4 
        notsharedN[6]= nn[6]
        notsharedN[1]= nn[1]
        notsharedN[2]= nn[2]
    elseif rel5
        notsharedN[1]= nn[1]
        notsharedN[2]= nn[2]
        notsharedN[3]= nn[3]
    elseif rel6
        notsharedN[2]= nn[2]
        notsharedN[3]= nn[3]
        notsharedN[4]= nn[4]
    end

    return notsharedN
end



























# Now I write a function which checks for the validity of a configuration, which also depends on the geometry of the configuration.
"""
    validConf2DGeneral(N,ind,edo,HPlist,dir,geometry)

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
                    end
                    
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
        
    
                    if dx ==1 && dy ==1 
                        ans=false
                        break
                    end
    
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
                    nnc=nearestNeighborsCoords2D(red,[x1,y1],triangular2D) # Nearest neigbors to our index `ind`.
                    
                    rel1= nnc[1] == periodicInd2D(red,[x2,y2])
                    rel2= nnc[2] == periodicInd2D(red,[x2,y2])
                    rel3= nnc[3] == periodicInd2D(red,[x2,y2])
                    rel4= nnc[4] == periodicInd2D(red,[x2,y2])
                    rel5= nnc[5] == periodicInd2D(red,[x2,y2])
                    rel6= nnc[6] == periodicInd2D(red,[x2,y2])
                    
                    overallcond= rel1 || rel2 || rel3 || rel4 || rel5 || rel6

                    if overallcond == false
                        ans=false
                        break
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
                    nnc=nearestNeighborsCoords2D(red,[x1,y1],triangular2D) # Nearest neigbors to our index `ind`.
                    
                    rel1= nnc[1] == periodicInd2D(red,[x2,y2])
                    rel2= nnc[2] == periodicInd2D(red,[x2,y2])
                    rel3= nnc[3] == periodicInd2D(red,[x2,y2])
                    rel4= nnc[4] == periodicInd2D(red,[x2,y2])
                    rel5= nnc[5] == periodicInd2D(red,[x2,y2])
                    rel6= nnc[6] == periodicInd2D(red,[x2,y2])
                    
                    overallcond= rel1 || rel2 || rel3 || rel4 || rel5 || rel6

                    if overallcond == false
                        ans=false
                        break
                    end
                end
            end
        end

    end

    return ans
end




























@enum NeighborsΔ::Int8 begin # Define nearest neighbors as enum type. The order mimicks the clockwise order.
        uprightΔ=1
        rightΔ=2
        downrightΔ=3
        downleftΔ=4
        leftΔ=5
        upleftΔ=6
    end 
"""
    countFirst2DΔ(N,edo,HPlist) 

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`,an array containing the sequence of H,P aminoacids `HPlist`, and
a geoetry; counts the number of possible pull moves for the first amino acid, recording the possible positions 
(for the triangle lattice).
"""
function countFirst2DΔ(N,edo,HPlist) 

    red=makeLattice(N,edo,HPlist)
    ind=1
    xplus,yplus=edo[ind+1,:]
    x,y=edo[ind,:]

    # Find out which spaces are free for `ind` to move into.
    freespaces=excludedNeighbors2D(red,[x,y],[xplus,yplus])

    nnc=nearestNeighborsCoords2D(red,[x,y],triangular2D)
    

    # Here I will store the possible coordinates for `ind`.
    coordinates1=[]
    for k in [uprightΔ,rightΔ,downrightΔ,downleftΔ,leftΔ,upleftΔ] 
        space=freespaces[Int(k)]
        if space == 0
            push!(coordinates1,nnc[Int(k)])
        end
    end

    numpull=length(coordinates1)
    # Now we have the possible coordinates for the first monomer, as well as the number of possible pull moves for the first 
    # amino acid. Fo easier use, I turn the arrays into matrices.

    coordinates1=transpose(hcat(coordinates1...))
    
    return (numpull,coordinates1)
end































"""
    countLast2DΔ(N,edo,HPlist)

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`,an array containing the sequence of H,P aminoacids `HPlist`, and
a geometry; counts the number of possible pull moves for the last amino acid, recording the possible positions 
(for the triangle lattice).
"""
function countLast2DΔ(N,edo,HPlist)
    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)
    ind=length(HPlist)
    xminus,yminus=edo[ind-1,:]
    x,y=edo[ind,:]

    # Find out which spaces are free for `ind` to move into.
    freespaces=excludedNeighbors2D(red,[x,y],[xminus,yminus])

    nnc=nearestNeighborsCoords2D(red,[x,y],triangular2D)
    

    # Here I will store the possible coordinates for `ind`.
    coordinates1=[]
    for k in [uprightΔ,rightΔ,downrightΔ,downleftΔ,leftΔ,upleftΔ]
        space=freespaces[Int(k)]
        if space == 0
            push!(coordinates1,nnc[Int(k)])
        end
    end

    numpull=length(coordinates1)
    # Now we have the possible coordinates for the first monomer, as well as the number of possible pull moves for the first 
    # amino acid. Fo easier use, I turn the arrays into matrices.

    coordinates1=transpose(hcat(coordinates1...))
    
    
    return (numpull,coordinates1)
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


    # Next we store the possible free spaces for `ind` when pulling backwards (as in moving only the previous monomers).
    freespaces1=sharedNeighbors2D(red,[x,y],[xplus,yplus])
    
    
    # Store possible free spaces when pulling forwards.
    freespaces2=sharedNeighbors2D(red,[x,y],[xminus,yminus])

    # Compute coordinates of nearest neighbors to `ind`.
    nnc=nearestNeighborsCoords2D(red,[x,y],triangular2D)
    

    # Now I write the possible coordinates for `ind`.
    coordinates1=[]
    for k in [uprightΔ,rightΔ,downrightΔ,downleftΔ,leftΔ,upleftΔ]
        space=freespaces1[Int(k)]
        if space == 0 
            push!(coordinates1,nnc[Int(k)])
        else
            continue
        end
    end

    coordinatesi=[]
    for k in [uprightΔ,rightΔ,downrightΔ,downleftΔ,leftΔ,upleftΔ]
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


























#Next, I write a function which counts all of the possible pull moves for a given state.
"""
    countpull2DΔ(N,edo,HPlist)

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; counts all of the possible pull moves. It also outputs the necessary
coordinates to perform one of the listed moves (for triangular geometry).
"""
function countpull2DΔ(N,edo,HPlist)

    # Count the moves and store the coordinates for the first monomer.
    t1,coords1=countFirst2DΔ(N,edo,HPlist)

    # Count the moves and store the coordinates for the last monomer.
    t2,coords3=countLast2DΔ(N,edo,HPlist)

    

    # Count the moves and store the coordinates for the middle amino acids. 
    t3=0
    tmid=0
    indexb=Int8[] # Stores the value of the middle indices that have non empty arrays of possible backwards pull moves.
    indexf=Int8[] # Stores the value of the middle indices that have non empty arrays of possible forwards pull moves.
    npullpindexb=Int8[] # Stores the number of backwards pull moves for each middle index.
    # Length may be shorter than the number of middle indices given that I only store the coordinates for non empty arrays.
    npullpindexf=Int8[] # Stores the number of forwards pull moves for each middle index.
    middlecoords1=[] # Contains the coordinates for `ind`, where `ind` is an index from the middle of the chain.
    middlecoords3=[] # Contains the coordinates for `ind` when the chain is pulled forwards.
    
    for j in 2:length(HPlist)-1
        tj,coords5,coords7=countMiddle2DΔ(N,j,edo,HPlist)
        # I have all the neccesary information. But first I need to check if the given index
        # posseses at least a possible pull move.
        t3=t3+tj
        if isempty(coords5) == false
            tmid=tmid+size(coords5)[1]
            push!(indexb,j)
            push!(npullpindexb,size(coords5)[1])
            push!(middlecoords1,coords5)
        end

        if isempty(coords7) == false
            push!(indexf,j)
            push!(npullpindexf,size(coords7)[1])
            push!(middlecoords3,coords7)
        end
    end

    totalpull=t1+t2+t3

    matrixb=hcat(indexb,npullpindexb)
    matrixf=hcat(indexf,npullpindexf)

    return (totalpull,tmid,matrixb,matrixf,coords1,coords3,middlecoords1,middlecoords3)
end
















# For convenience purposes, I write a function which counts the number of pull moves according to the geometry.
"""
    countpull2Dg(N,edo,HPlist,geometry)

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`, an array containing the sequence of H,P aminoacids 
`HPlist`, and a geometry; returns the number of possible pull moves.
"""
function totalpull2Dg(N,edo,HPlist,geometry)
    totalpull=0
    if geometry == square2D
        totalpull=countpull2DΔ(N,edo,HPlist)[1]
    elseif geometry == triangular2D
        totalpull=countpull2D(N,edo,HPlist)[1]
    end

    return totalpull
end


















"""
    pullMove2DΔ(N,edo,HPlist)  

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; chooses an index and performs a full pull move for the generated index 
(for the triangular geometry). 
"""
function pullMove2DΔ(N,edo,HPlist)  
    
    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)
    
    # `newedo` will contain the amino acids´ final positions.
    newedo=copy(edo)

    # Generate the list of possible pull moves.
    totalpull,tmid,matrixb,matrixf,coords1,coords3,middlecoords1,middlecoords3=countpull2DΔ(N,edo,HPlist)

    # Randomly choose one of the pull moves.
    m=rand(1:totalpull)
    
    # Make subdivisions according to the type of move.
    if isempty(coords1)
        s1=0
    else
        s1=size(coords1)[1] 
    end

    if isempty(coords3)
        s2=0
        s3=s1

    else
        s2=s1+1
        s3=s1+(size(coords3)[1])
    end

    s4=s3+1
    s5=s3+tmid
    s6=s5+1
    

    # Type of pull move depends on the value of `m`.
    if m ≤ s1 # First monomer is moved.
        ind=1
        singleMove2D(red,edo[ind,:],coords1[m,:]) # Move `ind`.
        newedo[ind,:]=coords1[m,:] # Record the change.

        # If all went well, a move has been done, now I have to check whether the current configuration is valid.
        stateconf=validConf2DGeneral(N,ind,newedo,HPlist,forwards,triangular2D)
        if stateconf == false
            for k in (ind+1):length(HPlist)
                if stateconf == false
                    op=edo[k,:] # Old position.
                    np=edo[k-1,:] # New position.
                    singleMove2D(red,op,np) 
                    newedo[k,:]=np 
                    stateconf=validConf2DGeneral(N,k,newedo,HPlist,forwards,triangular2D) # Check whether the new configuration is valid.
                else
                    break
                end
            end
        end





    elseif  s2 ≤ m ≤ s3
        mp=m-(s2-1) # `mp is the position` for the moves corresponding to the final monomer in the chain.
        ind=length(HPlist)
        singleMove2D(red,edo[ind,:],coords3[mp,:]) # Move `ind`.
        newedo[ind,:]=coords3[mp,:] # Record the change.

        # If all went well, a move has been done, now I have to check whether the current configuration is valid.
        stateconf=validConf2DGeneral(N,ind,newedo,HPlist,backwards,triangular2D)
        if stateconf == false
            for k in (ind-1):-1:1
                if stateconf == false
                    op=edo[k,:] # Old position.
                    np=edo[k+1,:] # New position.
                    singleMove2D(red,op,np) 
                    newedo[k,:]=np 
                    stateconf=validConf2DGeneral(N,k,newedo,HPlist,backwards,triangular2D) # Check whether the new configuration is valid.
                else
                    break
                end
            end
        end
        




    elseif s4 ≤ m ≤ s5
        mp=m-(s4-1)
        indm,ind,pos=middleInd2D(mp,matrixb)
        singleMove2D(red,edo[indm,:],middlecoords1[ind][pos,:])
        newedo[indm,:]=middlecoords1[ind][pos,:]

        # If all went well, a move has been done, now I have to check whether the current configuration is valid.
        stateconf=validConf2DGeneral(N,indm,newedo,HPlist,backwards,triangular2D)
        if stateconf == false
            for k in (indm-1):-1:1
                if stateconf == false
                    op=edo[k,:] # Old position.
                    np=edo[k+1,:] # New position.
                    singleMove2D(red,op,np) 
                    newedo[k,:]=np 
                    stateconf=validConf2DGeneral(N,k,newedo,HPlist,backwards,triangular2D) # Check whether the new configuration is valid.
                else
                    break
                end
            end
        end
        


    elseif m ≥ s6 
        mp=m-(s6-1)
        indm,ind,pos=middleInd2D(mp,matrixf)
        singleMove2D(red,edo[indm,:],middlecoords3[ind][pos,:])
        newedo[indm,:]=middlecoords3[ind][pos,:]

        # If all went well, a move has been performed, now I have to check whether the current configuration is valid.
        stateconf=validConf2DGeneral(N,indm,newedo,HPlist,forwards,triangular2D)
        if stateconf ==  false
            for k in (indm+1):length(HPlist)
                if stateconf == false
                    op=edo[k,:] # Old position.
                    np=edo[k-1,:] # New position.
                    singleMove2D(red,op,np) 
                    newedo[k,:]=np 
                    stateconf=validConf2DGeneral(N,k,newedo,HPlist,forwards,triangular2D) # Check whether the new configuration is valid.
                else
                    break
                end
            end
        end



    end
    
    newred=makeLattice(N,newedo,HPlist)
    
    return (newred,newedo,totalpull) # Returns everything neccesary to reproduce the final configuration and implement
    # the Metropolis algorithm.
end












"""
    pullMove2DGeneral(N,edo,HPlist,geometry)  

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`, an array containing 
the sequence of H,P aminoacids `HPlist`, and a geometry; chooses an index and performs a full pull move for 
the generated index.
"""
function pullMove2DGeneral(N,edo,HPlist,geometry) 
    if geometry == square2D
        newred,newedo,totalpull=pullMove2D(N,edo,HPlist)
    elseif geometry == triangular2D
        newred,newedo,totalpull=pullMove2DΔ(N,edo,HPlist)
    end

    return (newred,newedo,totalpull)
end












# infot=pullMove2DΔ(12,[[7 3];[6 3];[7 4];[8 5];[7 5];[7 6];[6 5];[5 5];[5 4];[6 4]],[1,-1,-1,1,-1,1,-1,1,-1,-1])


