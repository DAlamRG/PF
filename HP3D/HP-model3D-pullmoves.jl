# THis file contains the code necessary to perform a pull-move for the fcc geometry.



#Before doing anything, I define a structure which will encode the "protein" object.

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
    end








"""
    periodicInd3D(L,indices)

Given the size of a 3D array `L` and a couple of indices, returns the indices for the equivalent array with periodic boundary conditions.
"""
function periodicInd3D(L,indices)
    ix=iy=iz=L
    Ix=mod1(ix,lx)
    Iy=mod1(iy,ly)
    Iz=mod1(iz,lz)
        
    return [Ix,Iy,Iz]
end















"""
    periodicArr3D(L,indices)

    Given the size of a 3D array `L` and a couple of indices, returns the value of the equivalent array with periodic boundary conditions.
"""
function periodicArr3D(L,indices)
    ix=iy=iz=L
    Ix=mod1(ix,lx)
    Iy=mod1(iy,ly)
    Iz=mod1(iz,lz)
        
    return A[Ix,Iy,Iz]
end




















# Since I end up setting up the array which contains the protein sequence multiple times, I write a function to do 
#just that.
"""
    makeLattice3D(N,edo,HPlist)

Given a value for the lattice size `N`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; creates a 3D array containing the amino acid sequence.
"""
function makeLattice3D(N,edo,HPlist)
    red=zeros(Int8,(N,N,N))
    for k in 1:length(HPlist)
        x=edo[k,1]
        y=edo[k,2]
        z=edo[k,3]
        x,y,z=periodicInd3D(red,[x,y,z])
        red[x,y,z]=HPlist[k]
    end
    return red
end


















# Next, generalize the concept of nearest neighbors, according to the geometry.
"""
    nearestNeighbors3D(red,inds,geometry)

Given a 3D array `red`, and a couple of indices `inds`, and a geometry `geometry`; returns the values of the nearest neighbors to the given position. 
"""
function nearestNeighbors3D(red,inds,geometry)
    A=red
    x,y,z=inds
    x,y,z=periodicInd3D(red,[x,y,z])

    if geometry == cubic 
        nn=Int8[periodicArr3D(A,[x-1,y,z]),periodicArr3D(A,[x,y+1,z]),periodicArr3D(A,[x+1,y,z]),periodicArr3D(A,[x,y-1,z]),
        periodicArr3D(A,[x,y,z+1]),periodicArr3D(A,[x,y,z-1])] # Last two neighbors fall outside x-y plane.
    
    
    elseif geometry == hcp # First six neighbors are in the same x-y plane. The remaining six are in other planes.
        nn=Int8[periodicArr3D(A,[x-1,y,z]),periodicArr3D(A,[x,y+1,z]),periodicArr3D(A,[x+1,y+1,z]),periodicArr3D(A,[x+1,y,z])
        ,periodicArr3D(A,[x,y-1,z]),periodicArr3D(A,[x-1,y-1,z]),
        periodicArr3D(A,[x,y,z-1]),periodicArr3D(A,[x+1,y,z-1]),periodicArr3D(A,[x,y-1,z-1]),
        periodicArr3D(A,[x,y,z+1]),periodicArr3D(A,[x+1,y,z+1]),periodicArr3D(A,[x,y-1,z+1])]
    end
    
    return nn
end


















# Next, generalize the concept of nearest neighbors, according to the geometry.
"""
    nearestNeighborsCoords3D(red,inds,geometry)

Given a 3D array `red`, and a couple of indices `inds`, and a geometry `geometry`; returns the coordinates of the nearest neighbors 
to the given position. 
"""
function nearestNeighborsCoords3D(red,inds,geometry)
    A=red
    x,y,z=inds
    x,y,z=periodicInd3D(red,[x,y,z])

    if geometry == cubic 
        nnc=Vector{Int8}[periodicInd3D(A,[x-1,y,z]),periodicInd3D(A,[x,y+1,z]),periodicInd3D(A,[x+1,y,z]),periodicInd3D(A,[x,y-1,z]),
        periodicInd3D(A,[x,y,z+1]),periodicInd3D(A,[x,y,z-1])] # Last two neighbors fall outside x-y plane.
    
    
    elseif geometry == hcp # First six neighbors are in the same x-y plane. The remaining six are outside.
        nnc=Vector{Int8}[periodicInd3D(A,[x-1,y,z]),periodicInd3D(A,[x,y+1,z]),periodicInd3D(A,[x+1,y+1,z]),periodicInd3D(A,[x+1,y,z])
        ,periodicInd3D(A,[x,y-1,z]),periodicInd3D(A,[x-1,y-1,z]),
        periodicInd3D(A,[x,y,z-1]),periodicInd3D(A,[x+1,y,z-1]),periodicInd3D(A,[x,y-1,z-1]),
        periodicInd3D(A,[x,y,z+1]),periodicInd3D(A,[x+1,y,z+1]),periodicInd3D(A,[x,y-1,z+1])]
    end
    
    return nnc
end



















"""
    sharedNeighborsCoords3D(red,inds,indsp,geometry)

Given a 3D array `red`, a couple of indices `inds` and `indsp`, and a geometry; returns the coordinates for the empty shared 
neighbors by `inds,indsp` for the given 3D geometry. 
"""
function sharedNeighborsCoords3D(red,inds,indsp,geometry)
    nnc=nearestNeighborsCoords3D(red,inds,geometry) # Nearest neigbors to our index `ind`.
    nn=nearestNeighbors3D(red,inds,geometry)
    nncp=nearestNeighborsCoords3D(red,indsp,geometry)
    
    if geometry == cubic
        println("Hello, this function´s not complete yet")


    elseif geometry == hcp
        sharedNC=Vector{Int8}[] # Empty shared neighbor spaces will have a value of zero.
    
        for i in 1:12
            if (nnc[i] in nncp) && (nn[i] == 0)
                push!(sharedNC,nnc[i]) # These are the coordinates where we might move the index corresponding to `inds`
            end
        end
    end
    return sharedNC
end






















"""
    excludedNeighborsCoords3D(red,inds,indsp,geometry)

Given a 3D array `red`, a couple of indices `inds` and `indsp`, and a geometry; returns the coordinates for the empty not shared neighbor spaces by 
`inds,indsp` for the given 3D geometry.  
"""
function excludedNeighborsCoords3D(red,inds,indsp,geometry)
    nnc=nearestNeighborsCoords3D(red,inds,geometry) # Nearest neigbors to our index `ind`.
    nn=nearestNeighbors3D(red,inds,geometry)
    nncp=nearestNeighborsCoords3D(red,indsp,geometry)
    
    if geometry == cubic
        println("Hello, this function´s not complete yet")


    elseif geometry == hcp
        notsharedNC=Vector{Int8}[] # Empty shared neighbor spaces will have a value of zero.
        
        for i in 1:12
            if !(nnc[i] in nncp) && (nn[i] == 0) 
                push!(notsharedNC,nnc[i])
            end

        end

    end
    return notsharedNC
end





















# Now I write a function which checks for the validity of a configuration, which also depends on the geometry of the configuration.
"""
    validConf3D(N,ind,edo,HPlist,dir,geometry)

Given a 3D array of side `N`, a matrix encoding the aminoacids positions `edo`, an index `ind` , an array containing 
the sequence of H,P aminoacids `HPlist`, a direction `dir`, and a geometry; determines whether the protein structure is valid or not.
"""
function validConf3D(N,ind,edo,HPlist,dir,geometry)
    
    # I set up the protein within the array.
    red=makeLattice3D(N,edo,HPlist)

    # Next, I iterate over the protein´s vertices, computing the "distance" to the next vertex. The direction in which I check the
    # structure is determined by `dir`.
    ans=true

    if geometry ==  cubic # Check for validty in square geometry.
        println("Not yet finished")
    
    
    
    
    
    elseif geometry == hcp # Check the validty in triangular geometry.
        if dir == backwards
            if ind !=1
                for j in ind:-1:2
                    x1,y1,z1=edo[j,:]
                    x1,y1,z1=periodicInd3D(red,[x1,y1,z1]) # Make the indices periodic.
                    x2,y2,z2=edo[j-1,:]
                    x2,y2,z2=periodicInd3D(red,[x2,y2,z2]) # Make the indices periodic.
                    nnc=nearestNeighborsCoords3D(red,[x1,y1,z1],geometry) # Nearest neigbors to our index `ind`.
                    

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
                    x1,y1,z1=periodicInd3D(red,[x1,y1,z1]) # Make the indices periodic.
                    x2,y2,z2=edo[j+1,:]
                    x2,y2,z2=periodicInd3D(red,[x2,y2,z2]) # Make the indices periodic.
                    nnc=nearestNeighborsCoords3D(red,[x1,y1,z1],geometry) # Nearest neigbors to our index `ind`.
                    

                    if ([x2,y2,z2] in nnc) == false # Configuration is not valid.
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
    countFirst3D(N,edo,HPlist,geometry)

Given a 3D array size `N`, a matrix encoding the aminoacids positions `edo`,an array containing the sequence of H,P aminoacids `HPlist`, and
a geometry; counts the number of possible pull moves for the first amino acid, recording the possible positions.
"""
function countFirst3D(N,edo,HPlist,geometry)  
    
    if geometry == cubic
        println("Work in progress")
    
    elseif geometry == hcp
        red=makeLattice3D(N,edo,HPlist)
        ind=1
        xplus,yplus,zplus=edo[ind+1,:]
        x,y,z=edo[ind,:]

        # Find out which spaces are free for `ind` to move into.
        coordinates1=excludedNeighborsCoords3D(red,[x,y,z],[xplus,yplus,zplus],geometry)

        numpull=length(coordinates1)
        # Now we have the possible coordinates for the first monomer, as well as the number of possible pull moves for the first 
        # amino acid. Fo easier use, I turn the arrays into matrices.

        coordinates1=transpose(hcat(coordinates1...))

    end
    return (numpull,coordinates1)
end

















"""
    countLast3D(N,edo,HPlist,geometry) 

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`,an array containing the sequence of H,P aminoacids `HPlist`, and
a geometry; counts the number of possible pull moves for the first amino acid, recording the possible positions.
"""
function countLast3D(N,edo,HPlist,geometry) 

    red=makeLattice3D(N,edo,HPlist)
    ind=length(HPlist)
    xminus,yminus,zminus=edo[ind-1,:]
    x,y,z=edo[ind,:]

    if geometry == cubic
        println("Not yet done")

    elseif geometry == hcp
        # Find out which spaces are free for `ind` to move into.
        coordinates1=excludedNeighborsCoords3D(red,[x,y,z],[xminus,yminus,zminus],geometry)

        numpull=length(coordinates1)
        # Now we have the possible coordinates for the first monomer, as well as the number of possible pull moves for the first 
        # amino acid. Fo easier use, I turn the arrays into matrices.

        coordinates1=transpose(hcat(coordinates1...))

    end
        
    return (numpull,coordinates1)
end
















"""
    countMiddle3D(N,ind,edo,HPlist,geometry) 

Given a 3D array size `N`, an index `ind`, a matrix encoding the aminoacids positions `edo`, an array containing 
the sequence of H,P aminoacids `HPlist`, and a geometry; counts the number of possible pull moves for the middle amino acids, and 
stores the possible coordinates.
"""
function countMiddle3D(N,ind,edo,HPlist,geometry)
    
    # I set up the protein within the array.
    red=makeLattice3D(N,edo,HPlist)
    

    # Now we check wether a pull-move is possible for the selected vertex. A shared neighbor space between `ind` and either `ind+1` or `ind-1`
    # must be empty.
    # First I check where `ind+1,ind-1` are.
    xplus,yplus,zplus=edo[ind+1,:] # Position of `ind+1`
    xminus,yminus,zminus=edo[ind-1,:] # Position of `ind-1`
    x,y,z=edo[ind,:] # Position of `ind`


    # Next we store the possible free spaces for `ind` when pulling backwards (as in moving only the previous monomers).
    coordinates1=sharedNeighborsCoords3D(red,[x,y,z],[xplus,yplus,zplus],geometry)
    
    
    # Store possible free spaces when pulling forwards.
    coordinatesi=sharedNeighborsCoords3D(red,[x,y,z],[xminus,yminus,zminus],geometry)
    
    
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

















"""
    countpull3D(N,edo,HPlist,geometry)

Given a 3D array size `N`, a matrix encoding the aminoacids positions `edo`, an array containing 
the sequence of H,P aminoacids `HPlist`, and a geometry; counts all of the possible pull moves. It also outputs the necessary
coordinates to perform the listed moves.
"""
function countpull3D(N,edo,HPlist,geometry)

    # Count the moves and store the coordinates for the first monomer.
    t1,coords1=countFirst3D(N,edo,HPlist,geometry)

    # Count the moves and store the coordinates for the last monomer.
    t2,coords3=countLast3D(N,edo,HPlist,geometry)

    
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
        tj,coords5,coords7=countMiddle3D(N,j,edo,HPlist,geometry)
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



















"""
    singleMove(arr,pos1,pos2)

Given a 2D array `arr`, an initial position  `pos1` ([xᵢ,yᵢ]) and a final desired position `pos2` ([xₛ,yₛ]); moves the value of 
`arr[xᵢ,yᵢ]` to  `arr[xₛ,yₛ]`.
"""
function singleMove(arr,pos1,pos2)
    red=copy(arr)
    xi,yi,zi=pos1
    xi,yi,zi=periodicInd3D(red,[xi,yi,zi])
    xf,yf,zf=pos2
    xf,yf,zf=periodicInd3D(red,[xf,yf,zf])
    if pos1 != pos2
        red[xf,yf,zf]=arr[xi,yi,zf]
        red[xi,yi,zf]=0
    end
    return red
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
    pullMove3D(N,edo,HPlist,geometry)  

Given a 3D array size `N`, a matrix encoding the aminoacids positions `edo`, an array containing 
the sequence of H,P aminoacids `HPlist`, and a geometry; chooses an index and performs a full pull move for the generated index. 
"""
function pullMove3D(N,edo,HPlist,geometry)    
    
    # I set up the protein within the array.
    #red=makeLattice3D(N,edo,HPlist)
    
    # `newedo` will contain the amino acids´ final positions.
    newedo=copy(edo)

    # Generate the list of possible pull moves.
    totalpull,tmid,matrixb,matrixf,coords1,coords3,middlecoords1,middlecoords3=countpull3D(N,edo,HPlist,geometry)

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
        #singleMove(red,edo[ind,:],coords1[m,:]) # Move `ind`.
        newedo[ind,:]=coords1[m,:] # Record the change.

        # If all went well, a move has been done, now I have to check whether the current configuration is valid.
        stateconf=validConf3D(N,ind,newedo,HPlist,forwards,geometry)
        if stateconf == false
            for k in (ind+1):length(HPlist)
                if stateconf == false
                    #op=edo[k,:] # Old position.
                    np=edo[k-1,:] # New position.
                    #singleMove(red,op,np) 
                    newedo[k,:]=np 
                    stateconf=validConf3D(N,k,newedo,HPlist,forwards,geometry) # Check whether the new configuration is valid.
                else
                    break
                end
            end
        end





    elseif  s2 ≤ m ≤ s3
        mp=m-(s2-1) # `mp is the position` for the moves corresponding to the final monomer in the chain.
        ind=length(HPlist)
        #singleMove(red,edo[ind,:],coords3[mp,:]) # Move `ind`.
        newedo[ind,:]=coords3[mp,:] # Record the change.

        # If all went well, a move has been done, now I have to check whether the current configuration is valid.
        stateconf=validConf3D(N,ind,newedo,HPlist,backwards,geometry)
        if stateconf == false
            for k in (ind-1):-1:1
                if stateconf == false
                    #op=edo[k,:] # Old position.
                    np=edo[k+1,:] # New position.
                    #singleMove(red,op,np) 
                    newedo[k,:]=np 
                    stateconf=validConf3D(N,k,newedo,HPlist,backwards,geometry) # Check whether the new configuration is valid.
                else
                    break
                end
            end
        end
        




    elseif s4 ≤ m ≤ s5
        mp=m-(s4-1)
        indm,ind,pos=middleInd(mp,matrixb)
        #singleMove(red,edo[indm,:],middlecoords1[ind][pos,:])
        newedo[indm,:]=middlecoords1[ind][pos,:]

        # If all went well, a move has been done, now I have to check whether the current configuration is valid.
        stateconf=validConf3D(N,indm,newedo,HPlist,backwards,geometry)
        if stateconf == false
            for k in (indm-1):-1:1
                if stateconf == false
                    #op=edo[k,:] # Old position.
                    np=edo[k+1,:] # New position.
                    #singleMove(red,op,np) 
                    newedo[k,:]=np 
                    stateconf=validConf3D(N,k,newedo,HPlist,backwards,geometry) # Check whether the new configuration is valid.
                else
                    break
                end
            end
        end
        


    elseif m ≥ s6 
        mp=m-(s6-1)
        indm,ind,pos=middleInd(mp,matrixf)
        #singleMove(red,edo[indm,:],middlecoords3[ind][pos,:])
        newedo[indm,:]=middlecoords3[ind][pos,:]

        # If all went well, a move has been performed, now I have to check whether the current configuration is valid.
        stateconf=validConf3D(N,indm,newedo,HPlist,forwards,geometry)
        if stateconf ==  false
            for k in (indm+1):length(HPlist)
                if stateconf == false
                    #op=edo[k,:] # Old position.
                    np=edo[k-1,:] # New position.
                    #singleMove(red,op,np) 
                    newedo[k,:]=np 
                    stateconf=validConf3D(N,k,newedo,HPlist,forwards,geometry) # Check whether the new configuration is valid.
                else
                    break
                end
            end
        end



    end
    
    newred=makeLattice3D(N,newedo,HPlist)
    
    return (newred,newedo,totalpull) # Returns everything neccesary to reproduce the final configuration and implement
    # the Metropolis algorithm.
end
