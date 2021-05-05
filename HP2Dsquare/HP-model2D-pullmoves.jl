# Here I store all of the functions necessary to perform the so called pull-moves on a 2D array.

include("./visualizeHP2D.jl")



"""
    periodicInd2D(A,indices)

Given a 2D array `A` and a couple of indices, returns the indices for the equivalent array with periodic boundary conditions.
"""
function periodicInd2D(A,indices)
    lx,ly=size(A)
    ix,iy=indices
    Ix=mod1(ix,lx)
    Iy=mod1(iy,ly)
        
    return [Ix,Iy]
end

















"""
    periodicArr2D(A,indices)

Given a 2D array `A` and a couple of indices, returns the value of the equivalent array with periodic boundary conditions.
"""
function periodicArr2D(A,indices)
    lx,ly=size(A)
    ix,iy=indices
    Ix=mod1(ix,lx)
    Iy=mod1(iy,ly)
        
    return A[Ix,Iy]
end













# Since I end up setting up the array which contains the protein sequence multiple times, I write a function to do 
#just that.
"""
    makeLattice(N,edo,HPlist)

Given a value for the lattice size `N`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; creates a 2D array containing the amino acid sequence.
"""
function makeLattice(N,edo,HPlist)
    red=zeros(Int8,(N,N))
    for k in 1:length(HPlist)
        x=edo[k,1]
        y=edo[k,2]
        x,y=periodicInd2D(red,[x,y])
        red[x,y]=HPlist[k]
    end
    return red
end
























# Now, the task is to write a function which determines wether a given protein sequence is valid.
"""
    validConf(N,ind,edo,HPlist)

Given a 2D array of side `N`, a matrix encoding the aminoacids positions `edo`, an index `ind` , and an array containing 
the sequence of H,P aminoacids `HPlist`; determines whether the protein structure is valid or not.
"""
function validConf(N,ind,edo,HPlist)
    
    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)
    
    # Next, I iterate over the protein´s vertices, computing the "distance" to the next vertix.
    ans=true
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
            elseif dy > 1
                if ycond == false
                    ans=false
                    break
                end
            end
        end
    end
return ans
end
            


















"""
    singleMove2D(arr,pos1,pos2)

Given a 2D array `arr`, an initial position  `pos1` ([xᵢ,yᵢ]) and a final desired position `pos2` ([xₛ,yₛ]); moves the value of 
`arr[xᵢ,yᵢ]` to  `arr[xₛ,yₛ]`.
"""
function singleMove2D(arr,pos1,pos2)
    red=copy(arr)
    xi,yi=pos1
    xi,yi=periodicInd2D(red,[xi,yi])
    xf,yf=pos2
    xf,yf=periodicInd2D(red,[xf,yf])
    if pos1 != pos2
        red[xf,yf]=arr[xi,yi]
        red[xi,yi]=0
    end
    return red
end




















"""
    function makecross2D(red,inds)

Given a 2D array `red`, and a couple of indices `inds`; returns the values of the nearest neighbors to the given position,
in clockwise order. 
"""
function makecross2D(red,inds)
    A=red
    x,y=inds
    x,y=periodicInd2D(red,[x,y])

    cross=Int8[periodicArr2D(A,[x-1,y]),periodicArr2D(A,[x,y+1]),periodicArr2D(A,[x+1,y]),periodicArr2D(A,[x,y-1])]
    return cross
end
    

















"""
    distance2D(red,inds1,inds2)

 Given an array `red`, two pairs of indices `inds1,inds2`; return the distance between the given positions.
"""
function distance2D(red,inds1,inds2)
    x1,y1=inds1
    x2,y2=inds2
    x1,y1=periodicInd2D(red,inds1)
    x2,y2=periodicInd2D(red,inds2)

    N=size(red)[1]

    xcond1= x1 == 1 && x2 == N
    xcond2= x1 == N && x2 == 1
    xcond= xcond1 || xcond2 
            
    ycond1= y1 == 1 && y2 == N
    ycond2= y1 == N && y2 == 1
    ycond= ycond1 || ycond2

    if xcond
        dx=1
    else
        dx=abs(x1-x2)
    end

    if ycond
        dy=1
    else
        dy=abs(y1-y2)
    end

    d=dx+dy
    return (d,dx,dy)
end

















"""
    validConfEnd2D(red,ind,inds,edo,HPlist)

Given a 2D array `red`, a position `ind`, a couple of indices `inds`, a matrix encoding the aminoacids positions `edo`,
 and an array containing the sequence of H,P aminoacids `HPlist`; checks whether the position `inds` for the selected amino acid
 is valid. The necessary condition is for the end amino acid not to be adjacent to the chain.
"""
function validConfEnd2D(red,ind,inds,edo,HPlist)
    
    x,y=inds
    x,y=periodicInd2D(red,[x,y])

    # `ind` can only be `1` or `length(HPlist)`, meaning this function only checks for the validity of a configuration for the ends of the chain.

    value = true
    if ind == 1
        for k in 3:length(HPlist) # Iterate over the mebers of the chain not directly linked to `ind`.
            xk,yk=edo[k,:]
            xk,yk=periodicInd2D(red,[xk,yk])
            dist=distance2D(red,[x,y],[xk,yk])[1]
            if dist == 1 # End monomer is adjacent to the chain, configuration is not valid.
                value=false
                break
            end
        end
    elseif ind == length(HPlist)
        for k in length(HPlist)-2:-1:1
            xk,yk=edo[k,:]
            xk,yk=periodicInd2D(red,[xk,yk])
            dist=distance2D(red,[x,y],[xk,yk])[1]
            if dist == 1 
                value=false
                break
            end

    end
    return value
end



















"""
    countFirst2D(N,edo,HPlist)

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; counts the number of possible pull moves for the first amino acid, 
recording the possible positions.
"""
function countFirst2D(N,edo,HPlist)
    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)
    ind=1
    xplus,yplus=edo[ind+1,:]
    x,y=edo[ind,:]
    dx,dy=distance2D(red,[xplus,yplus],[x,y])[2:3]
    rel1= dx > 0 && dy == 0
    rel2= dx == 0 && dy < 0
    rel3= dx < 0 && dy == 0
    rel4= dx == 0 && dy > 0 # Conditions single out possible adjacent spaces to move the second amino acid to.

    @enum Neighbors::Int8 begin # Define nearest neighbors as enum type.
        up=1
        right=2
        down=3
        left=4
    end

    vs=makecross2D(red,[x,y]) # These are the nearest neighbors to `ind`.
    vscoords=[[x-1,y],[x,y+1],[x+1,y],[x,y-1]]

    crossv1=makecross2D(red,[x-1,y])
    crossv2=makecross2D(red,[x,y+1])
    crossv3=makecross2D(red,[x+1,y])
    crossv4=makecross2D(red,[x,y-1]) # Nearest neighbors of the nearest neighbors.

    indices1=Int8[] # Here I will store the coordinates for `ind+1`.
    indices2=Int8[] # Here I store the possible coordinates for `ind`.
    if rel1 
        if vs[2] == 0 # `ind+1` might move into `vs[2]`.
            for k in [right,down]
                if crossv2[Int(k)] == 0
                    push!(indices1,2)
                    push!(indices2,k)
                end
            end
        elseif vs[3] == 0
            for k in [right,down,left]
                if crossv3[Int(k)] == 0
                    push!(indices1,3)
                    push!(indices2,k)
                end
            end
        elseif vs[4] == 0
            for k in [down,left]
                if crossv4[Int(k)] == 0
                    push!(indices1,4)
                    push!(indices2,k)
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
        elseif vs[3] == 0
            for k in [down,left]
                if crossv3[Int(k)] == 0
                    push!(indices1,3)
                    push!(indices2,k)
                end
            end
        elseif vs[4] == 0
            for k in [up,down,left]
                if crossv4[Int(k)] == 0
                    push!(indices1,4)
                    push!(indices2,k)
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
        elseif vs[2] == 0
            for k in [up,right]
                if crossv2[Int(k)] == 0
                    push!(indices1,2)
                    push!(indices2,k)
                end
            end
        elseif vs[4] == 0
            for k in [up,left]
                if crossv4[Int(k)] == 0
                    push!(indices1,4)
                    push!(indices2,k)
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
        elseif vs[2] == 0
            for k in [up,right,down]
                if crossv2[Int(k)] == 0
                    push!(indices1,2)
                    push!(indices2,k)
                end
            end
        elseif vs[3] == 0
            for k in [right,down]
                if crossv3[Int(k)] == 0
                    push!(indices1,3)
                    push!(indices2,k)
                end
            end

    end
    # I now have two arrays with the information of where it is possible to move the first amino acids. First, I convert the numbers or enum types to coordinates.
    # Then, making use of `validConfEnd2D` I make sure that the end amino acid `ind` doesn´t end up in the vicinity of another 
    # amino acid besides `ind+1`.

    coordinates1=[]
    coordinates2=[]

    for el in indices1
        if el == 1
            push!(coordinates1,vscoords[el])
        elseif el == 2
            push!(coordinates1,vscoords[el])
        elseif el == 3
            push!(coordinates1,vscoords[el])
        elseif el == 4
            push!(coordinates1,vscoords[el])
        end
    end

    for k in 1:length(indices2)
        el=indices2[k]
        xv,yv=coordinates1[k]
        if el == up
            push!(coordinates2,[xv-1,yv])
        elseif el == right
            push!(coordinates2,[xv,yv+1])
        elseif el == down
            push!(coordinates2,[xv+1,yv])
        elseif el == left
            push!(coordinates2,[xv,yv-1])
        end
    end

    # Now I have the coordinates for `ind` and `ind+1`. I proceed to check wheter the new configuration would be valid.

    for k in 1:length(indices1)
        xa,ya=coordinates2[k] # Possible new coordinates for the first amino acid.

        if validConfEnd2D(red,1,[xa,ya],edo,HPlist) != true
            deleteat!(coordinates1,k)
            deleteat!(coordinates2,k) # Delete possible coordinates, doesn´t satisfy the condition.
        end
    end

    numpull=length(coordinates1)
    # Now we have the possible coordinates for the first two monomers, as well as the number of possible pull moves for the first 
    # amino acid. Fo easier use, i turn the arrays `coordinates1,coordinates2` into a `numpull×2` matrices.

    coordinates1=transpose(hcat(coordinates1...))
    coordinates2=transpose(hcat(coordinates2...))
    return (numpull,coordinates1,coordinates2)
end

























"""
    countLast2D(N,edo,HPlist)

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; counts the number of possible pull moves for the last amino acid, 
recording the possible positions.
"""
function countLast2D(N,edo,HPlist)
    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)
    ind=length(HPlist)
    xminus,yminus=edo[ind-1,:]
    x,y=edo[ind,:]
    dx,dy=distance2D(red,[xminus,yminus],[x,y])[2:3]
    rel1= dx > 0 && dy == 0
    rel2= dx == 0 && dy < 0
    rel3= dx < 0 && dy == 0
    rel4= dx == 0 && dy > 0 # Conditions single out possible adjacent spaces to move the second to last amnino acid to.

    @enum Neighbors::Int8 begin # Define nearest neighbors as enum type.
        up=1
        right=2
        down=3
        left=4
    end

    vs=makecross2D(red,[x,y]) # These are the nearest neighbors to `ind`.
    vscoords=[[x-1,y],[x,y+1],[x+1,y],[x,y-1]]

    crossv1=makecross2D(red,[x-1,y])
    crossv2=makecross2D(red,[x,y+1])
    crossv3=makecross2D(red,[x+1,y])
    crossv4=makecross2D(red,[x,y-1]) # Nearest neighbors of the nearest neighbors.

    indices1=Int8[] # Here I will store the coordinates for `ind-1`.
    indices2=Int8[] # Here I store the possible coordinates for `ind`.
    if rel1 
        if vs[2] == 0 # `ind+1` might move into `vs[2]`.
            for k in [right,down]
                if crossv2[Int(k)] == 0
                    push!(indices1,2)
                    push!(indices2,k)
                end
            end
        elseif vs[3] == 0
            for k in [right,down,left]
                if crossv3[Int(k)] == 0
                    push!(indices1,3)
                    push!(indices2,k)
                end
            end
        elseif vs[4] == 0
            for k in [down,left]
                if crossv4[Int(k)] == 0
                    push!(indices1,4)
                    push!(indices2,k)
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
        elseif vs[3] == 0
            for k in [down,left]
                if crossv3[Int(k)] == 0
                    push!(indices1,3)
                    push!(indices2,k)
                end
            end
        elseif vs[4] == 0
            for k in [up,down,left]
                if crossv4[Int(k)] == 0
                    push!(indices1,4)
                    push!(indices2,k)
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
        elseif vs[2] == 0
            for k in [up,right]
                if crossv2[Int(k)] == 0
                    push!(indices1,2)
                    push!(indices2,k)
                end
            end
        elseif vs[4] == 0
            for k in [up,left]
                if crossv4[Int(k)] == 0
                    push!(indices1,4)
                    push!(indices2,k)
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
        elseif vs[2] == 0
            for k in [up,right,down]
                if crossv2[Int(k)] == 0
                    push!(indices1,2)
                    push!(indices2,k)
                end
            end
        elseif vs[3] == 0
            for k in [right,down]
                if crossv3[Int(k)] == 0
                    push!(indices1,3)
                    push!(indices2,k)
                end
            end

    end
    # I now have two arrays with the information of where it is possible to move the last amino acids. First, I convert the numbers 
    # or enum types to coordinates. Then, making use of `validConfEnd2D` I make sure that the end amino acid `ind` doesn´t 
    # end up in the vicinity of another amino acid besides `ind-1`.

    coordinates1=[]
    coordinates2=[]

    for el in indices1
        if el == 1
            push!(coordinates1,vscoords[el])
        elseif el == 2
            push!(coordinates1,vscoords[el])
        elseif el == 3
            push!(coordinates1,vscoords[el])
        elseif el == 4
            push!(coordinates1,vscoords[el])
        end
    end

    for k in 1:length(indices2)
        el=indices2[k]
        xv,yv=coordinates1[k]
        if el == up
            push!(coordinates2,[xv-1,yv])
        elseif el == right
            push!(coordinates2,[xv,yv+1])
        elseif el == down
            push!(coordinates2,[xv+1,yv])
        elseif el == left
            push!(coordinates2,[xv,yv-1])
        end
    end

    # Now I have the coordinates for `ind` and `ind-1`. I proceed to check wheter the new configuration would be valid.

    for k in 1:length(indices1)
        xa,ya=coordinates2[k] # Possible new coordinates for the first amino acid.

        if validConfEnd2D(red,length(HPlist),[xa,ya],edo,HPlist) != true
            deleteat!(coordinates1,k)
            deleteat!(coordinates2,k) # Delete possible coordinates, they don´t satisfy the condition.
        end
    end

    numpull=length(coordinates1)
    # Now we have the possible coordinates for the last two monomers, as well as the number of possible pull moves for the last 
    # amino acid. Fo easier use, i turn the arrays `coordinates1,coordinates2` into a `numpull×2` matrices.

    coordinates1=transpose(hcat(coordinates1...))
    coordinates2=transpose(hcat(coordinates2...))
    return (numpull,coordinates1,coordinates2)
end
































"""
    countMiddle2D(N,ind,edo,HPlist)

Given a 2D array size `N`, an index `ind`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; counts the number of possible pull moves for the middle amino acids, and stores the possible coordinates. 
"""
function countMiddle2D(N,ind,edo,HPlist)
    
    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)
    display(red)
    
    # `newedo` will contain the amino acids´ final positions.
    newedo=copy(edo)
        
    # Now we check wether a pull-move is possible for the selected vertex. A diagonal space must be empty. Also, the empty
    # diagonal space must be adjacent to the amino acid correpsonding to `ind+1`.
    # First I check where `ind+1` is.
    xplus=edo[ind+1,1]
    yplus=edo[ind+1,2] # Position of `ind+1`
    x=edo[ind,1]
    y=edo[ind,2] # Position of `ind`
    dx=x-xplus
    dy=y-yplus
    rel1= dx > 0 && dy == 0
    rel2= dx == 0 && dy < 0
    rel3= dx < 0 && dy == 0
    rel4= dx == 0 && dy > 0 # Conditions single out possible diagonal spaces to move to.

    # Next we store the possible diagonal spaces.
    diagspaces=ones(Int8,4)
    if rel1 
        diagspaces[1]=red[x-1,y+1]
        diagspaces[4]=red[x-1,y-1]
    elseif rel2 
        diagspaces[1]=red[x-1,y+1]
        diagspaces[2]=red[x+1,y+1]
    elseif rel3 
        diagspaces[2]=red[x+1,y+1]
        diagspaces[3]=red[x+1,y-1]
    elseif rel4 
        diagspaces[3]=red[x+1,y-1]
        diagspaces[4]=red[x-1,y-1]
    end


    v1=red[x-1,y]
    v2=red[x,y+1]
    v3=red[x+1,y]
    v4=red[x,y-1]  # Possible sites to move `ind-1` to.

    cond1= edo[ind-1,:] == [x-1,y]
    cond2= edo[ind-1,:] == [x,y+1]
    cond3= edo[ind-1,:] == [x+1,y]
    cond4= edo[ind-1,:] == [x,y-1] # `edo[ind-1,:]` is the current location on the array `red` of `ind-1`.
    
    condi1= (v1 == 0) || cond1 # Either the place is empty, or `ind-1` is already there.
    condi2= (v2 == 0) || cond2
    condi3= (v3 == 0) || cond3
    condi4= (v4 == 0) || cond4
    condis=Bool[condi1 || condi2, condi2 || condi3, condi3 || condi4, condi4 || condi1]



    indices=Int8[]
    for k in 1:4
        space=diagspaces[k]
        # If the diagonal space is empty, and the adjacent space is either empty or occupied by `ind-1` ,
        # pull move might be possible.
        if (space == 0 && condis[k] ) 
            push!(indices,k)
        else
            continue
        end
    end


    # I have an array with the possible spaces to move to. Now I have to check whether there is space for the amino acid 
    # correpsonding to `ind-1` to move to.
    if isempty(indices) == false


        coordinates1=[] # Possible coordinates for `ind`
        coordinates2=[] # Possible coordinates for `ind-1`.
        xys=zeros(Int8,2) # To be filled with the position `ind-1` might be needed to be moved to.

        for el in indices

            if el == 1
                if condi1  # Either `v1` is empty or it is already occupied by `ind-1`. A pull move is all but assured.
                    singleMove2D(red,[x,y],[x-1,y+1]) # Make the move.
                    newedo[ind,:]=[x-1,y+1] # Record the change.
                    xys[:]=[x-1,y]
                    break
                elseif condi2 
                    singleMove2D(red,[x,y],[x-1,y+1])
                    newedo[ind,:]=[x-1,y+1]
                    xys[:]=[x,y+1]
                    break
                end


            elseif el == 2
                if condi2 
                    singleMove2D(red,[x,y],[x+1,y+1])
                    newedo[ind,:]=[x+1,y+1]
                    xys[:]=[x,y+1]
                    break
                elseif condi3 
                    singleMove2D(red,[x,y],[x+1,y+1])
                    newedo[ind,:]=[x+1,y+1]
                    xys[:]=[x+1,y]
                    break  
                end


            elseif el == 3 
                if condi3 
                    singleMove2D(red,[x,y],[x+1,y-1])
                    newedo[ind,:]=[x+1,y-1]
                    xys[:]=[x+1,y]
                    break
                elseif condi4 
                    singleMove2D(red,[x,y],[x+1,y-1])
                    newedo[ind,:]=[x+1,y-1]
                    xys[:]=[x,y-1]
                    break
                end


            elseif el == 4 
                if condi4 
                    singleMove2D(red,[x,y],[x-1,y-1])
                    newedo[ind,:]=[x-1,y-1]
                    xys[:]=[x,y-1]
                    break
                elseif condi1 
                    singleMove2D(red,[x,y],[x-1,y-1])
                    newedo[ind,:]=[x-1,y-1]
                    ys[:]=[x-1,y]
                    break
                end
            end
        end
        
        # Now I have in `xys` the potential new position for `ind-1`. If the current configuration is valid, I don´t make the move
        # and the full pull move has been succesfully completed.
        if validConf(N,ind,newedo,HPlist) != true
            singleMove2D(red,edo[ind-1,:],xys) # Makes the move.
            newedo[ind-1,:]=xys[:] # Records the change.
        end

    end
    

    # If all went well, a move has been done, now I have to check whether the current configuration is valid.
    stateconf=validConf(N,ind-1,newedo,HPlist) 
    for k in ind:-1:3
        if stateconf == false
            op=edo[k-2,:] # Old position.
            np=edo[k,:] # New position.
            singleMove2D(red,op,np) 
            newedo[k-2,:]=np 
            stateconf=validConf(N,k-2,newedo,HPlist) # Check whether the new configuration is valid.
        else
            break
        end
    end

    return (red,newedo) # Returns everything neccesary to reproduce the final configuration.
end



















"""
    pullMove2D(N,ind,edo,HPlist)

Given a 2D array size `N`, an index `ind`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; counts the number of possible pull moves for the middle amino acids, and stores the possible coordinates. 
"""
function pullMove2D(N,ind,edo,HPlist)
    
    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)
    display(red)
    
    # `newedo` will contain the amino acids´ final positions.
    newedo=copy(edo)
        
    # Now we check wether a pull-move is possible for the selected vertex. A diagonal space must be empty. Also, the empty
    # diagonal space must be adjacent to the amino acid correpsonding to `ind+1`.
    # First I check where `ind+1` is.
    xplus=edo[ind+1,1]
    yplus=edo[ind+1,2] # Position of `ind+1`
    x=edo[ind,1]
    y=edo[ind,2] # Position of `ind`
    dx=x-xplus
    dy=y-yplus
    rel1= dx > 0 && dy == 0
    rel2= dx == 0 && dy < 0
    rel3= dx < 0 && dy == 0
    rel4= dx == 0 && dy > 0 # Conditions single out possible diagonal spaces to move to.

    # Next we store the possible diagonal spaces.
    diagspaces=ones(Int8,4)
    if rel1 
        diagspaces[1]=red[x-1,y+1]
        diagspaces[4]=red[x-1,y-1]
    elseif rel2 
        diagspaces[1]=red[x-1,y+1]
        diagspaces[2]=red[x+1,y+1]
    elseif rel3 
        diagspaces[2]=red[x+1,y+1]
        diagspaces[3]=red[x+1,y-1]
    elseif rel4 
        diagspaces[3]=red[x+1,y-1]
        diagspaces[4]=red[x-1,y-1]
    end


    v1=red[x-1,y]
    v2=red[x,y+1]
    v3=red[x+1,y]
    v4=red[x,y-1]  # Possible sites to move `ind-1` to.

    cond1= edo[ind-1,:] == [x-1,y]
    cond2= edo[ind-1,:] == [x,y+1]
    cond3= edo[ind-1,:] == [x+1,y]
    cond4= edo[ind-1,:] == [x,y-1] # `edo[ind-1,:]` is the current location on the array `red` of `ind-1`.
    
    condi1= (v1 == 0) || cond1 # Either the place is empty, or `ind-1` is already there.
    condi2= (v2 == 0) || cond2
    condi3= (v3 == 0) || cond3
    condi4= (v4 == 0) || cond4
    condis=Bool[condi1 || condi2, condi2 || condi3, condi3 || condi4, condi4 || condi1]



    indices=Int8[]
    for k in 1:4
        space=diagspaces[k]
        # If the diagonal space is empty, and the adjacent space is either empty or occupied by `ind-1` ,
        # pull move might be possible.
        if (space == 0 && condis[k] ) 
            push!(indices,k)
        else
            continue
        end
    end


    # I have an array with the possible spaces to move to. Now I have to check whether there is space for the amino acid 
    # correpsonding to `ind-1` to move to.
    if isempty(indices) == false

        xys=zeros(Int8,2) # To be filled with the position `ind-1` might be needed to be moved to.

        for el in indices

            if el == 1
                if condi1  # Either `v1` is empty or it is already occupied by `ind-1`. A pull move is all but assured.
                    singleMove2D(red,[x,y],[x-1,y+1]) # Make the move.
                    newedo[ind,:]=[x-1,y+1] # Record the change.
                    xys[:]=[x-1,y]
                    break
                elseif condi2 
                    singleMove2D(red,[x,y],[x-1,y+1])
                    newedo[ind,:]=[x-1,y+1]
                    xys[:]=[x,y+1]
                    break
                end


            elseif el == 2
                if condi2 
                    singleMove2D(red,[x,y],[x+1,y+1])
                    newedo[ind,:]=[x+1,y+1]
                    xys[:]=[x,y+1]
                    break
                elseif condi3 
                    singleMove2D(red,[x,y],[x+1,y+1])
                    newedo[ind,:]=[x+1,y+1]
                    xys[:]=[x+1,y]
                    break  
                end


            elseif el == 3 
                if condi3 
                    singleMove2D(red,[x,y],[x+1,y-1])
                    newedo[ind,:]=[x+1,y-1]
                    xys[:]=[x+1,y]
                    break
                elseif condi4 
                    singleMove2D(red,[x,y],[x+1,y-1])
                    newedo[ind,:]=[x+1,y-1]
                    xys[:]=[x,y-1]
                    break
                end


            elseif el == 4 
                if condi4 
                    singleMove2D(red,[x,y],[x-1,y-1])
                    newedo[ind,:]=[x-1,y-1]
                    xys[:]=[x,y-1]
                    break
                elseif condi1 
                    singleMove2D(red,[x,y],[x-1,y-1])
                    newedo[ind,:]=[x-1,y-1]
                    ys[:]=[x-1,y]
                    break
                end
            end
        end
        
        # Now I have in `xys` the potential new position for `ind-1`. If the current configuration is valid, I don´t make the move
        # and the full pull move has been succesfully completed.
        if validConf(N,ind,newedo,HPlist) != true
            singleMove2D(red,edo[ind-1,:],xys) # Makes the move.
            newedo[ind-1,:]=xys[:] # Records the change.
        end

    end
    

    # If all went well, a move has been done, now I have to check whether the current configuration is valid.
    stateconf=validConf(N,ind-1,newedo,HPlist) 
    for k in ind:-1:3
        if stateconf == false
            op=edo[k-2,:] # Old position.
            np=edo[k,:] # New position.
            singleMove2D(red,op,np) 
            newedo[k-2,:]=np 
            stateconf=validConf(N,k-2,newedo,HPlist) # Check whether the new configuration is valid.
        else
            break
        end
    end

    return (red,newedo) # Returns everything neccesary to reproduce the final configuration.
end







