# Generalized pull-moves for both dimensiona and the three geometries.

# First, declare the relevant structures for both dimensions.


@enum geometries::Int8 begin # First, define the type of geometries as an enum type.
    square2D=1
    triangular2D=2
    cubic=3
    fcc=4
end


struct Protein
    edo :: Matrix{Int64} # This encodes the amino acidds positions within the array.
    HPlist :: Array{Int8,1} # This encodes the protein´s amino acid sequence.
    geometry :: Enum # This defines the type of geometry that the protein is embedded in.
end


@enum directions::Int8 begin
        forwards=1
        backwards=2
        nonetaken=3
    end








"""
    periodicInd(A,indices,dim)

Given a 3D array `A`, a couple/triad of indices and a dimension `dim`; returns the indices for the equivalent array 
with periodic boundary conditions.
"""
function periodicInd(A,indices,dim)
    if dim == 2
        lx,ly=size(A)
        ix,iy=indices
        Ix=mod1(ix,lx)
        Iy=mod1(iy,ly)
        return [Ix,Iy]
    else
        lx,ly,lz=size(A)
        ix,iy,iz=indices
        Ix=mod1(ix,lx)
        Iy=mod1(iy,ly)
        Iz=mod1(iz,lz)
        return [Ix,Iy,Iz]
    end
end











"""
    periodicArr(A,indices,dim)

Given a 3D array `A`, a couple/triad of indices and a dimension `dim`; returns the value of the equivalent array 
with periodic boundary conditions.
"""
function periodicArr(A,indices,dim)
    if dim == 2
        lx,ly=size(A)
        ix,iy=indices
        Ix=mod1(ix,lx)
        Iy=mod1(iy,ly)
        return A[Ix,Iy]
    else
        lx,ly,lz=size(A)
        ix,iy,iz=indices
        Ix=mod1(ix,lx)
        Iy=mod1(iy,ly)
        Iz=mod1(iz,lz)
        return A[Ix,Iy,Iz]
    end
end










# Since I end up setting up the array which contains the protein sequence multiple times, I write a function to do 
#just that.
"""
    makeLattice(N,edo,HPlist)

Given a value for the lattice size `N`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; creates a 2D/3D array containing the amino acid sequence.
"""
function makeLattice(N,edo,HPlist)
    if length(edo[k,:]) == 3
        red=zeros(Int8,(N,N,N))
        for k in 1:length(HPlist)
            x=edo[k,1]
            y=edo[k,2]
            z=edo[k,3]
            x,y,z=periodicInd(red,[x,y,z],3)
            red[x,y,z]=HPlist[k]
        end    
        return red
    else
        red=zeros(Int8,(N,N))
        for k in 1:length(HPlist)
            x=edo[k,1]
            y=edo[k,2]
            x,y=periodicInd(red,[x,y],2)
            red[x,y]=HPlist[k]
        end    
        return red
    end
end














# Next, generalize the concept of nearest neighbors, according to the geometry.
"""
    nearestNeighbors(red,inds,geometry)

Given a 2D/3D array `red`, and a couple of indices `inds`, and a geometry `geometry`; returns the values of the nearest 
neighbors to the given position. 
"""
function nearestNeighbors(red,inds,geometry)
    A=red

    if geometry == cubic 
        x,y,z=inds
        x,y,z=periodicInd(red,[x,y,z],3)
        nn=Int8[periodicArr(A,[x-1,y,z],3),periodicArr(A,[x,y+1,z],3),periodicArr(A,[x+1,y,z],3),periodicArr(A,[x,y-1,z],3),
        periodicArr(A,[x,y,z+1],3),periodicArr(A,[x,y,z-1],3)] # Last two neighbors fall outside x-y plane.
        return nn
    

    elseif geometry == fcc # Fcc geometry has 12 topological nearest neighbors.
        i,j,k=inds
        i,j,k=periodicInd(red,[i,j,k],3)
        nn=Int8[periodicArr(A,[i+1,j,k],3),periodicArr(A,[i-1,j,k],3)
        ,periodicArr(A,[i,j+1,k],3),periodicArr(A,[i,j-1,k],3)
        ,periodicArr(A,[i,j,k+1],3),periodicArr(A,[i,j,k-1],3)
        ,periodicArr(A,[i+1,j-1,k],3),periodicArr(A,[i-1,j+1,k],3)
        ,periodicArr(A,[i,j+1,k-1],3),periodicArr(A,[i,j-1,k+1],3)
        ,periodicArr(A,[i-1,j,k+1],3),periodicArr(A,[i+1,j,k-1],3)]
        return nn
    

    elseif geometry == square2D
        x,y=inds
        x,y=periodicInd(red,[x,y],2)
        nn=Int8[periodicArr(A,[x-1,y],2),periodicArr(A,[x,y+1],2),periodicArr(A,[x+1,y],2),periodicArr(A,[x,y-1],2)]
        return nn


    elseif geometry == triangular2D
        x,y=inds
        x,y=periodicInd(red,[x,y],2)
        nn=Int8[periodicArr(A,[x-1,y],2),periodicArr(A,[x,y+1],2)
        ,periodicArr(A,[x+1,y+1],2),periodicArr(A,[x+1,y],2)
        ,periodicArr(A,[x,y-1],2),periodicArr(A,[x-1,y-1],2)]
        return nn
    end

end
















# Next, generalize the concept of nearest neighbors, according to the geometry.
"""
    nearestNeighborsCoords(red,inds,geometry)

Given a 2D/3D array `red`, and a couple of indices `inds`, and a geometry `geometry`; returns the coordinates of the nearest neighbors 
to the given position. 
"""
function nearestNeighborsCoords(red,inds,geometry)
    A=red

    if geometry == cubic 
        x,y,z=inds
        x,y,z=periodicInd(red,[x,y,z],3)
        nnc=Vector{Int8}[periodicInd(A,[x-1,y,z],3),periodicInd(A,[x,y+1,z],3),periodicInd(A,[x+1,y,z],3),periodicInd(A,[x,y-1,z],3),
        periodicInd(A,[x,y,z+1],3),periodicInd(A,[x,y,z-1],3)] # Last two neighbors fall outside x-y plane.
        return nnc
    
    
    elseif geometry == fcc # First six neighbors are in the same x-y plane. The remaining six are outside.
        i,j,k=inds
        i,j,k=periodicInd(red,[i,j,k],3)
        nnc=Vector{Int8}[periodicInd(A,[i+1,j,k],3),periodicInd(A,[i-1,j,k],3)
        ,periodicInd(A,[i,j+1,k],3),periodicInd(A,[i,j-1,k],3)
        ,periodicInd(A,[i,j,k+1],3),periodicInd(A,[i,j,k-1],3)
        ,periodicInd(A,[i+1,j-1,k],3),periodicInd(A,[i-1,j+1,k],3)
        ,periodicInd(A,[i,j+1,k-1],3),periodicInd(A,[i,j-1,k+1],3)
        ,periodicInd(A,[i-1,j,k+1],3),periodicInd(A,[i+1,j,k-1],3)]
        return nnc
    

    elseif geometry == square2D
        x,y=inds
        x,y=periodicInd(red,[x,y],2)
        nnc=Vector{Int8}[periodicInd(A,[x-1,y],2),periodicInd(A,[x,y+1],2),periodicInd(A,[x+1,y],2),periodicInd(A,[x,y-1],2)]
        return nnc


    elseif geometry == triangular2D
        x,y=inds
        x,y=periodicInd(red,[x,y],2)
        nnc=Vector{Int8}[periodicInd(A,[x-1,y],2),periodicInd(A,[x,y+1],2)
        ,periodicInd(A,[x+1,y+1],2),periodicInd(A,[x+1,y],2)
        ,periodicInd(A,[x,y-1],2),periodicInd(A,[x-1,y-1],2)]
        return nnc
    end
end

















"""
    sharedNeighborsCoords(red,inds,indsp,geometry)

Given a 2D/3D array `red`, a couple of indices `inds` and `indsp`, and a geometry; returns the coordinates for the empty shared topological
neighbors by `inds,indsp` for the given geometry. 
"""
function sharedNeighborsCoords(red,inds,indsp,geometry)
    nnc=nearestNeighborsCoords(red,inds,geometry) # Coordiantes of nearest neigbors to our index `ind`.
    nn=nearestNeighbors(red,inds,geometry) # Value of nearest neighbors to our index.
    nncp=nearestNeighborsCoords(red,indsp,geometry) # Coordiantes of nearest neigbors to our index `indsp`.
    
    if geometry == cubic
        neighbornumber = 6
    elseif geometry == fcc
        neighbornumber = 12
    elseif geometry == square2D
        neighbornumber = 4
    elseif geometry == triangular2D
        neighbornumber = 6
    end

    sharedNC=Vector{Int8}[] # Empty shared neighbor spaces will have a value of zero.
    for i in 1:neighbornumber
        if (nnc[i] in nncp) && (nn[i] == 0)
            push!(sharedNC,nnc[i]) # These are the coordinates where we might move the index corresponding to `inds`
        end
    end

    return sharedNC
end
















"""
    excludedNeighborsCoords(red,inds,indsp,geometry)

Given a 2D/3D array `red`, a couple of indices `inds` and `indsp`, and a geometry; returns the coordinates for the empty not-shared neighbor spaces by 
`inds,indsp` for the given 3D geometry.  
"""
function excludedNeighborsCoords(red,inds,indsp,geometry)
    nnc=nearestNeighborsCoords(red,inds,geometry) # Coordiantes of nearest neigbors to our index `ind`.
    nn=nearestNeighbors(red,inds,geometry) # Value of nearest neighbors to our index.
    nncp=nearestNeighborsCoords(red,indsp,geometry) # Coordiantes of nearest neigbors to our index `indsp`.

    if geometry == cubic
        neighbornumber = 6
    elseif geometry == fcc
        neighbornumber = 12
    elseif geometry == square2D
        neighbornumber = 4
    elseif geometry == triangular2D
        neighbornumber = 6
    end
    

    notsharedNC=Vector{Int8}[] # Empty not shared neighbor spaces will have a value of zero.
    for i in 1:neighbornumber
        if !(nnc[i] in nncp) && (nn[i] == 0)
            push!(notsharedNC,nnc[i]) # These are the coordinates where we might move the index corresponding to `inds`
        end
    end

    return notsharedNC
end















# Now I write a function which checks for the validity of a configuration, which also depends on the geometry of the configuration.
"""
    validConf(N,ind,edo,HPlist,dir,geometry)

Given a 2D/3D array of side `N`, a matrix encoding the aminoacids positions `edo`, an index `ind` , an array containing 
the sequence of H,P aminoacids `HPlist`, a direction `dir`, and a geometry; determines whether the protein structure is valid or not.
"""
function validConf(N,ind,edo,HPlist,dir,geometry)

    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)

    # Next, I iterate over the protein´s vertices, cheacking wether each pair of positions is in each other´s list of nearest neighbors.
    # The direction in which I check the structure is determined by `dir`.
    ans=true

    
    if (geometry == fcc) || (geometry == cubic) # Check the validty in fcc geometry.
        if dir == backwards
            if ind !=1
                for j in ind:-1:2
                    x1,y1,z1=edo[j,:]
                    x1,y1,z1=periodicInd(red,[x1,y1,z1],3) # Make the indices periodic.
                    x2,y2,z2=edo[j-1,:]
                    x2,y2,z2=periodicInd(red,[x2,y2,z2],3) # Make the indices periodic.
                    nnc=nearestNeighborsCoords(red,[x1,y1,z1],geometry) # Nearest neigbors to the position `[x1,y1,z1]`.
                    

                    if ([x2,y2,z2] in nnc) == false # Configuration is not valid.
                        ans=false
                        break
                    end
                end
            end
    
    
        elseif dir == forwards
            if ind != length(HPlist)
                for j in ind:length(HPlist)-1
                    x1,y1,z1=edo[j,:]
                    x1,y1,z1=periodicInd(red,[x1,y1,z1],3) # Make the indices periodic.
                    x2,y2,z2=edo[j-1,:]
                    x2,y2,z2=periodicInd(red,[x2,y2,z2],3) # Make the indices periodic.
                    nnc=nearestNeighborsCoords(red,[x1,y1,z1],geometry) # Nearest neigbors to the position `[x1,y1,z1]`.
                    

                    if ([x2,y2,z2] in nnc) == false # Configuration is not valid.
                        ans=false
                        break
                    end
                end
            end

        end

    elseif (geometry == square2D) || (geometry == triangular2D) # Check the validty in square2D geometry.
        if dir == backwards
            if ind !=1
                for j in ind:-1:2
                    x1,y1=edo[j,:]
                    x1,y1=periodicInd(red,[x1,y1],2) # Make the indices periodic.
                    x2,y2=edo[j-1,:]
                    x2,y2=periodicInd(red,[x2,y2],2) # Make the indices periodic.
                    nnc=nearestNeighborsCoords(red,[x1,y1],geometry) # Nearest neigbors to the position `[x1,y1]`.
                    

                    if ([x2,y2] in nnc) == false # Configuration is not valid.
                        ans=false
                        break
                    end
                end
            end
    
    
        elseif dir == forwards
            if ind != length(HPlist)
                for j in ind:length(HPlist)-1
                    x1,y1=edo[j,:]
                    x1,y1=periodicInd(red,[x1,y1],2) # Make the indices periodic.
                    x2,y2=edo[j-1,:]
                    x2,y2=periodicInd(red,[x2,y2],2) # Make the indices periodic.
                    nnc=nearestNeighborsCoords(red,[x1,y1],geometry) # Nearest neigbors to the position `[x1,y1]`.
                    

                    if ([x2,y2] in nnc) == false # Configuration is not valid.
                        ans=false
                        break
                    end
                end
            end

        end

    end


    return ans
end















"""
    countFirst(N,edo,HPlist,geometry)

Given a 2D/3D array size `N`, a matrix encoding the aminoacids positions `edo`,an array containing the sequence of H,P aminoacids `HPlist`, and
a geometry; counts the number of possible pull moves for the first amino acid, recording the possible positions.
"""
function countFirst(N,edo,HPlist,geometry)  
    
    if geometry == cubic
        println("Work in progress")

    elseif geometry == square2D

    
    elseif (geometry == fcc) || (geometry == triangular2D)
        red=makeLattice(N,edo,HPlist)
        ind=1

        # Find out which spaces are free for `ind` to move into.
        coordinates1=excludedNeighborsCoords(red,edo[ind,:],edo[ind+1,:],geometry)

        numpull=length(coordinates1) # Compute the number of pull-moves.
        # Now we have the possible coordinates for the first monomer, as well as the number of possible pull moves for the first 
        # amino acid. Fo easier use, I turn the arrays into matrices.

        coordinates1=transpose(hcat(coordinates1...))
        return (numpull,coordinates1)
    end
end















"""
    countLast(N,edo,HPlist,geometry) 

Given a 2D/3D array size `N`, a matrix encoding the aminoacids positions `edo`,an array containing the sequence of H,P aminoacids `HPlist`, and
a geometry; counts the number of possible pull moves for the first amino acid, recording the possible positions.
"""
function countLast(N,edo,HPlist,geometry) 

    if geometry == cubic
        println("Work in progress")

    elseif geometry == square2D

    
    elseif (geometry == fcc) || (geometry == triangular2D)
        red=makeLattice(N,edo,HPlist)
        ind=length(HPlist)

        # Find out which spaces are free for `ind` to move into.
        coordinates1=excludedNeighborsCoords(red,edo[ind,:],edo[ind-1,:],geometry)

        numpull=length(coordinates1) # Compute the number of pull-moves.
        # Now we have the possible coordinates for the first monomer, as well as the number of possible pull moves for the first 
        # amino acid. Fo easier use, I turn the arrays into matrices.

        coordinates1=transpose(hcat(coordinates1...))
        return (numpull,coordinates1)
    end
end
















"""
    countMiddle(N,ind,edo,HPlist,geometry) 

Given a 2D/3D array size `N`, an index `ind`, a matrix encoding the aminoacids positions `edo`, an array containing 
the sequence of H,P aminoacids `HPlist`, and a geometry; counts the number of possible pull moves for the middle amino acids, and 
stores the possible coordinates.
"""
function countMiddle(N,ind,edo,HPlist,geometry)
    
    if (geometry == fcc) || (geometry == triangular2D)
        # I set up the protein within the array.
        red=makeLattice(N,edo,HPlist)
    
        # Now we check wether a pull-move is possible for the selected vertex. A shared neighbor space between `ind` and either `ind+1` or `ind-1`
        # must be empty.

        # Next we store the possible free spaces for `ind` when pulling backwards (as in moving only the previous monomers).
        coordinates1=sharedNeighborsCoords3D(red,edo[ind,:],edo[ind+1,:],geometry)
    
        # Store possible free spaces when pulling forwards.
        coordinatesi=sharedNeighborsCoords3D(red,edo[ind,:],edo[ind-1,:],geometry)
    
    
        # Now I have in `coordinates1,coordinatesi` the potential new positions for `ind`. The length of these arrays is 
        # the number of possible pull moves por the given index.
        numpull1=length(coordinates1)
        numpull2=length(coordinatesi)
        numpull=numpull1+numpull2

        # I turn the 1 dimensional arrays into matrices for easier use.
        coordinates1=transpose(hcat(coordinates1...))
        coordinatesi=transpose(hcat(coordinatesi...))

        return (numpull,coordinates1,coordinatesi) # Returns everything neccesary to reproduce the final configuration.
    elseif geometry == square2D
        # Work in progress.
    end
end
















"""
    countpull(N,edo,HPlist,geometry)

Given a 2D/3D array size `N`, a matrix encoding the aminoacids positions `edo`, an array containing 
the sequence of H,P aminoacids `HPlist`, and a geometry; counts all of the possible pull moves. It also outputs the necessary
coordinates to perform the listed moves.
"""
function countpull(N,edo,HPlist,geometry)

    if (geometry == fcc) || (geometry == triangular2D)
        # Count the moves and store the coordinates for the first monomer.
        t1,coords1=countFirst(N,edo,HPlist,geometry)

        # Count the moves and store the coordinates for the last monomer.
        t2,coords3=countLast(N,edo,HPlist,geometry)

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
            tj,coords5,coords7=countMiddle(N,j,edo,HPlist,geometry)
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

    elseif  geometry == square2D
        # Work in progress.
    end

end














# Before writing the function that performs the actual move, I need a function which returns the index being pullled
# in one of the `middlecoords`matrices from my function `countpull3D`.
"""
    middleInd(mp,matrix)

Given a number `mp`, and a matrix `matrix` containing the values of the middle indices with possible pull moves in 
the first column and the number of pull moves for said indices in the second column;returns the corresponding 
index `ind` being moved and the position in one of my `middlecoords[ind][position,:]` matrices.
"""
function middleInd(mp,matrix)
    indexd=matrix[:,1] # This array contains the value of the indices with viable pull moves.
    npullpindex=matrix[:,2] # This array contains the number of pull moves for each index in `indexd`.
    l=length(npullpindex)
    nb=ones(Int8,l) # `nb` will store the number of accumulated pull moves for each of the indices.
    nb[1]=npullpindex[1]
    for i in 2:l # Fill `nb`.
        nb[i]=npullpindex[i]+nb[i-1]
    end

    indm=0 # This is the index being pulled.
    ind=0 # This is the index for the `middlecoords[ind][pos,:]` matrices.
    pos=0 # This is the position of the right coordinates.
    
    
    for i in 1:l
        if mp ≤ nb[i]
            indm=indexd[i]
            ind=i
            break
        end
    end

    if mp ≤ nb[1]
        pos=mp
    else
        pos=mp-nb[ind-1]
    end
    
    return (indm,ind,pos)
end


















"""
    pullMove(N,edo,HPlist,geometry)  

Given a 2D/3D array size `N`, a matrix encoding the aminoacids positions `edo`, an array containing 
the sequence of H,P aminoacids `HPlist`, and a geometry; chooses an index and performs a full pull move for the generated index. 
The function also returns the final state `newedo`, a number of possible pull moves `totalpull` for the 
initial configuration, the position of the pulled monomer `indpulled`, and the new position for the pulled index `newcoord`.
"""
function pullMove(N,edo,HPlist,geometry)   


    if geometry == fcc || geometry == triangular2D

        # `newedo` will contain the amino acids´ final positions.
        newedo=copy(edo)

        # Generate the list of possible pull moves.
        totalpull,tmid,matrixb,matrixf,coords1,coords3,middlecoords1,middlecoords3=countpull(N,edo,HPlist,geometry)

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
            newedo[ind,:]=coords1[m,:] # Record the change.
            indpulled=ind # Function returns the index pulled in the move.
            newcoord=newedo[ind,:] # Function returns the new coordinate for the pulled index.
            dirpulled= forwards
    
            # If all went well, a move has been done, now I have to check whether the current configuration is valid.
            stateconf=validConf(N,ind,newedo,HPlist,forwards,geometry)
            if stateconf == false
                for k in (ind+1):length(HPlist)
                    if stateconf == false
                        np=edo[k-1,:] # New position.
                        newedo[k,:]=np 
                        stateconf=validConf(N,k,newedo,HPlist,forwards,geometry) # Check whether the new configuration is valid.
                    else
                        break
                    end
                end
            end
    
    
    
    
        elseif  s2 ≤ m ≤ s3
            mp=m-(s2-1) # `mp is the position` for the moves corresponding to the final monomer in the chain.
            ind=length(HPlist)
            newedo[ind,:]=coords3[mp,:] # Record the change.
            indpulled=ind
            newcoord=newedo[ind,:]
            dirpulled= backwards
    
            # If all went well, a move has been done, now I have to check whether the current configuration is valid.
            stateconf=validConf(N,ind,newedo,HPlist,backwards,geometry)
            if stateconf == false
                for k in (ind-1):-1:1
                    if stateconf == false
                        np=edo[k+1,:] # New position.
                        newedo[k,:]=np 
                        stateconf=validConf(N,k,newedo,HPlist,backwards,geometry) # Check whether the new configuration is valid.
                    else
                        break
                    end
                end
            end
            
    
    
    
        elseif s4 ≤ m ≤ s5
            mp=m-(s4-1)
            indm,ind,pos=middleInd(mp,matrixb)
            newedo[indm,:]=middlecoords1[ind][pos,:]
            indpulled=indm
            newcoord=newedo[indm,:]
            dirpulled= backwards
    
            # If all went well, a move has been done, now I have to check whether the current configuration is valid.
            stateconf=validConf(N,indm,newedo,HPlist,backwards,geometry)
            if stateconf == false
                for k in (indm-1):-1:1
                    if stateconf == false
                        np=edo[k+1,:] # New position.
                        newedo[k,:]=np 
                        stateconf=validConf(N,k,newedo,HPlist,backwards,geometry) # Check whether the new configuration is valid.
                    else
                        break
                    end
                end
            end
            
    
    
        elseif m ≥ s6 
            mp=m-(s6-1)
            indm,ind,pos=middleInd(mp,matrixf)
            newedo[indm,:]=middlecoords3[ind][pos,:]
            indpulled=indm
            newcoord=newedo[indm,:]
            dirpulled= forwards
    
            # If all went well, a move has been performed, now I have to check whether the current configuration is valid.
            stateconf=validConf(N,indm,newedo,HPlist,forwards,geometry)
            if stateconf ==  false
                for k in (indm+1):length(HPlist)
                    if stateconf == false
                        np=edo[k-1,:] # New position. 
                        newedo[k,:]=np 
                        stateconf=validConf(N,k,newedo,HPlist,forwards,geometry) # Check whether the new configuration is valid.
                    else
                        break
                    end
                end
            end
    
    
    
        end
    
        return (newedo,totalpull,indpulled,newcoord,dirpulled) # Returns everything neccesary to reproduce the final configuration and implement
        # the Metropolis algorithm.


    elseif geometry == square2D
        # Work in progress

    end
        
    
end



















"""
    reconstructStates(N,edo,HPlist,pulledindices,dirs,newcoords,geometry)
    
Given a lattice size `N`, a matrix encoding the aminoacids initial position `edo`, an aminoacid sequence `HPlist` ,an array containing the 
pulled indices `pulledindices`, an array containing the direction in which the chain was pulled `dirs`, a matrix containing the new positions 
for the pulled indices `newcoords`, a geometry; returns a multidimensional array containing the states at each point of the Metropolis scheme.
"""
function reconstructStates(N,edo,HPlist,pulledindices,dirs,newcoords,geometry)

    if geometry == fcc || geometry == triangular2D

        if geometry == fcc
            dim = 3
        else
            dim = 2
        end
        reconstructedSates=zeros(Int8,(length(HPlist),dim,length(pulledindices)+1))
        reconstructedSates[:,:,1]=edo
        for l in 2:length(pulledindices)+1
            ind=pulledindices[l-1]
    
            if ind == 0 # Pull move did not happen.
                reconstructedSates[:,:,l]=copy(reconstructedSates[:,:,l-1]) # State is unchanged
    
            else
                newedo=copy(reconstructedSates[:,:,l-1])
                refedo=copy(reconstructedSates[:,:,l-1])
                newedo[ind,:]=newcoords[l-1,:]
    
                if ind == 1
                    stateconf=validConf(N,ind,newedo,HPlist,forwards,geometry)
                    if stateconf == false
                        for k in (ind+1):length(HPlist)
                            if stateconf == false
                                np=refedo[k-1,:] # New position.
                                newedo[k,:]=np 
                                stateconf=validConf(N,k,newedo,HPlist,forwards,geometry) # Check whether the new configuration is valid.
                            else
                                break
                            end
                        end
                    end
                    reconstructedSates[:,:,l]=newedo
    
    
    
                elseif ind == length(HPlist)
                    stateconf=validConf(N,ind,newedo,HPlist,backwards,geometry)
                    if stateconf == false
                        for k in (ind-1):-1:1
                            if stateconf == false
                                np=refedo[k+1,:] # New position.
                                newedo[k,:]=np 
                                stateconf=validConf(N,k,newedo,HPlist,backwards,geometry) # Check whether the new configuration is valid.
                            else
                                break
                            end
                        end
                    end
                    reconstructedSates[:,:,l]=newedo
    
    
    
                else
                    dir=dirs[l-1]
    
                    if dir == backwards
                        stateconf=validConf(N,ind,newedo,HPlist,dir,geometry)
                        if stateconf == false
                            for k in (ind-1):-1:1
                                if stateconf == false
                                    np=refedo[k+1,:] # New position.
                                    newedo[k,:]=np 
                                    stateconf=validConf(N,k,newedo,HPlist,dir,geometry) # Check whether the new configuration is valid.
                                else
                                    break
                                end
                            end
                        end
                        reconstructedSates[:,:,l]=newedo
    
                    elseif dir == forwards
                        stateconf=validConf(N,ind,newedo,HPlist,dir,geometry)
                        if stateconf ==  false
                            for k in (ind+1):length(HPlist)
                                if stateconf == false
                                    np=refedo[k-1,:] # New position.
                                    newedo[k,:]=np 
                                    stateconf=validConf(N,k,newedo,HPlist,dir,geometry) # Check whether the new configuration is valid.
                                else
                                    break
                                end
                            end
                        end
                        reconstructedSates[:,:,l]=newedo
                    end
    
                    
                end
    
            end
    
        end

        return reconstructedSates


    elseif geometry == square2D
        # work in progress

    end
end


