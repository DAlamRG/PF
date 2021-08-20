# Here I store all of the functions necessary to perform the so called pull-moves on a 2D array.








#Before doing anything, I define a structure which will encode the "protein" object.

@enum geometries::Int8 begin # First, define the type of geometries as an enum type.
    square2D=1
    triangular2D=2
end



struct Protein2D
    edo :: Matrix{Int64} # This encodes the amino acidds positions within the array.
    HPlist :: Array{Int8,1} # This encodes the protein´s amino acid sequence.
    geometry :: Enum # This defines the type of geometry that the protein is embedded in.
end












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




















@enum directions::Int8 begin
        forwards=1
        backwards=2
        nonetaken=3
    end
# Now, the task is to write a function which determines whether a given protein sequence is valid. This function works 
# only for the middle indices.
"""
    validConf(N,ind,edo,HPlist,dir)

Given a 2D array of side `N`, a matrix encoding the aminoacids positions `edo`, an index `ind` , an array containing 
the sequence of H,P aminoacids `HPlist`, and a direction `dir`; determines whether the protein structure is valid or not.
"""
function validConf(N,ind,edo,HPlist,dir)
    
    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)

    # Next, I iterate over the protein´s vertices, computing the "distance" to the next vertex. The direction in which I check the
    # structure is determined by `dir`.
    ans=true

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
    makecross2D(red,inds)

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





















# I believe this function is not neccesary. I leave it here just in case it comes in handy later.
"""
    validConfEnd2D(red,ind,inds,edo,HPlist)

Given a 2D array `red`, a position `ind`, a couple of indices `inds`, a matrix encoding the aminoacids positions `edo`,
and an array containing the sequence of H,P aminoacids `HPlist`; checks whether the position `inds` for the selected amino acid
is valid. The necessary condition is for the end amino acid not to be adjacent to the chain.
"""
function validConfEnd2D(red,ind,inds,edo,HPlist)
    x,y=inds
    x,y=periodicInd2D(red,[x,y])

    # `ind` can only be `1` or `length(HPlist)`, meaning this function only checks for the validity of a configuration 
    # at the ends of the chain.

    value = true
    if ind == 1
        for k in 3:length(HPlist) # Iterate over the elements of the chain not directly linked to `ind`.
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

    end
    return value
end
























@enum Neighbors::Int8 begin # Define nearest neighbors as enum type.
        up=1
        right=2
        down=3
        left=4
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
    
    rel1= periodicInd2D(red,[x-1,y]) == periodicInd2D(red,[xminus,yminus])
    rel2= periodicInd2D(red,[x,y+1]) == periodicInd2D(red,[xminus,yminus])
    rel3= periodicInd2D(red,[x+1,y]) == periodicInd2D(red,[xminus,yminus])
    rel4= periodicInd2D(red,[x,y-1]) == periodicInd2D(red,[xminus,yminus]) # Conditions single out possible adjacent spaces to move the second to last amnino acid to.

    vs=makecross2D(red,[x,y]) # These are the nearest neighbors to `ind`.
    vscoords=[periodicInd2D(red,[x-1,y]),periodicInd2D(red,[x,y+1]),periodicInd2D(red,[x+1,y]),periodicInd2D(red,[x,y-1])] # coordinates
    # of nearest neighbors

    crossv1=makecross2D(red,[x-1,y])
    crossv2=makecross2D(red,[x,y+1])
    crossv3=makecross2D(red,[x+1,y])
    crossv4=makecross2D(red,[x,y-1]) # Nearest neighbors of the nearest neighbors.

    indices1=Int8[] # Here I will store the coordinates for `ind-1`.
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
    # I now have two arrays with the information of where it is possible to move the last amino acids. First, I convert the numbers 
    # or enum types to coordinates. Then, making use of `validConfEnd2D` I make sure that the end amino acid `ind` doesn´t 
    # end up in the vicinity of another amino acid besides `ind-1`.

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
    # Now we have the possible coordinates for the last two monomers, as well as the number of possible pull moves for the last 
    # amino acid. Fo easier use, I turn the arrays `coordinates1,coordinates2` into a `numpull×2` matrices.

    coordinates1=transpose(hcat(coordinates1...))
    coordinates2=transpose(hcat(coordinates2...))
    return (numpull,coordinates1,coordinates2)
end
























"""
    diagSpaces2D(red,inds,indsp)

Given a 2D array `red`, a couple of indices `inds` and `indsp`; returns the possible diagonal spaces to which `ind` might be moved to. 
"""
function diagSpaces2D(red,inds,indsp)
    diagspaces=ones(Int8,4) # Possible diagonal spaces will have a value of zero.
    x,y=periodicInd2D(red,inds) # Position of amino acid we wanr to move.
    
    
    rel1= periodicInd2D(red,[x-1,y]) == periodicInd2D(red,indsp)
    rel2= periodicInd2D(red,[x,y+1]) == periodicInd2D(red,indsp)
    rel3= periodicInd2D(red,[x+1,y]) == periodicInd2D(red,indsp)
    rel4= periodicInd2D(red,[x,y-1]) == periodicInd2D(red,indsp) # Possible relative positions between `inds,indsp`.

    if rel1 # Fill `diagspaces` according to the relative positions.
        diagspaces[1]= periodicArr2D(red,[x-1,y+1])
        diagspaces[4]= periodicArr2D(red,[x-1,y-1])
    elseif rel2 
        diagspaces[1]= periodicArr2D(red,[x-1,y+1])
        diagspaces[2]= periodicArr2D(red,[x+1,y+1])
    elseif rel3 
        diagspaces[2]= periodicArr2D(red,[x+1,y+1])
        diagspaces[3]= periodicArr2D(red,[x+1,y-1])
    elseif rel4 
        diagspaces[3]= periodicArr2D(red,[x+1,y-1])
        diagspaces[4]= periodicArr2D(red,[x-1,y-1])
    end

    return diagspaces
end




















"""
    countMiddle2D(N,ind,edo,HPlist)

Given a 2D array size `N`, an index `ind`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; counts the number of possible pull moves for the middle amino acids, and 
stores the possible coordinates. 
"""
function countMiddle2D(N,ind,edo,HPlist)
    
    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)
    

    # Now we check wether a pull-move is possible for the selected vertex. A diagonal space must be empty. Also, the empty
    # diagonal space must be adjacent to the amino acid correpsonding to `ind+1` or to `ind-1`.
    # First I check where `ind+1,ind-1` are.
    xplus,yplus=edo[ind+1,:] # Position of `ind+1`
    xminus,yminus=edo[ind-1,:] # Position of `ind-1`
    x,y=edo[ind,:] # Position of `ind`


    # Next we store the possible diagonal spaces for `ind` when we wish to pull backwards.
    diagspaces1=diagSpaces2D(red,[x,y],[xplus,yplus])
    
    
    # Store possible diagonal spaces when we wish to pull forwards.
    diagspaces2=diagSpaces2D(red,[x,y],[xminus,yminus])
    


    v1=periodicArr2D(red,[x-1,y])
    v2=periodicArr2D(red,[x,y+1])
    v3=periodicArr2D(red,[x+1,y])
    v4=periodicArr2D(red,[x,y-1])  # Possible sites to move `ind-1,ind+1` to.

    cond1= periodicInd2D(red,[xminus,yminus]) == periodicInd2D(red,[x-1,y])
    cond2= periodicInd2D(red,[xminus,yminus]) == periodicInd2D(red,[x,y+1])
    cond3= periodicInd2D(red,[xminus,yminus]) == periodicInd2D(red,[x+1,y])
    cond4= periodicInd2D(red,[xminus,yminus]) == periodicInd2D(red,[x,y-1]) # `edo[ind-1,:]` is the current location on the array `red` of `ind-1`.

    condit1= periodicInd2D(red,[xplus,yplus]) == periodicInd2D(red,[x-1,y])
    condit2= periodicInd2D(red,[xplus,yplus]) == periodicInd2D(red,[x,y+1])
    condit3= periodicInd2D(red,[xplus,yplus]) == periodicInd2D(red,[x+1,y])
    condit4= periodicInd2D(red,[xplus,yplus]) == periodicInd2D(red,[x,y-1]) # `edo[ind+1,:]` is the current location on the array `red` of `ind+1`.
    
    condi1= (v1 == 0) || cond1 # Either the place is empty, or `ind-1` is already there.
    condi2= (v2 == 0) || cond2
    condi3= (v3 == 0) || cond3
    condi4= (v4 == 0) || cond4
    condis=Bool[condi1 || condi2, condi2 || condi3, condi3 || condi4, condi4 || condi1]

    conditio1= (v1 == 0) || condit1 # Either the place is empty, or `ind+1` is already there.
    conditio2= (v2 == 0) || condit2
    conditio3= (v3 == 0) || condit3
    conditio4= (v4 == 0) || condit4
    condists=Bool[conditio1 || conditio2, conditio2 || conditio3, conditio3 || conditio4, conditio4 || conditio1]


    indices=Int8[]
    for k in 1:4
        space=diagspaces1[k]
        # If the diagonal space is empty, and the adjacent space is either empty or occupied by `ind-1` ,
        # pull move might be possible.
        if (space == 0 && condis[k] ) 
            push!(indices,k)
        else
            continue
        end
    end

    indices2=Int8[]
    for k in 1:4
        space=diagspaces2[k]
        # If the diagonal space is empty, and the adjacent space is either empty or occupied by `ind+1` ,
        # pull move might be possible.
        if (space == 0 && condists[k] ) 
            push!(indices2,k)
        else
            continue
        end
    end
    
    



    # I have an array with the possible spaces to move to. Now I have to check whether there is space for the amino acid 
    # correpsonding to `ind-1` to move to.
    coordinates1=[] # Possible coordinates for `ind`
    coordinates2=[] # Possible coordinates for `ind-1`.

    diag1=periodicInd2D(red,[x-1,y+1])
    diag2=periodicInd2D(red,[x+1,y+1])
    diag3=periodicInd2D(red,[x+1,y-1])
    diag4=periodicInd2D(red,[x-1,y-1]) # These are the coordinates I am going to push into `coordinates1`.

    # The next array contains the coordinates I will push into `coordinates2`.
    vscoords=[periodicInd2D(red,[x-1,y]),periodicInd2D(red,[x,y+1]),periodicInd2D(red,[x+1,y]),periodicInd2D(red,[x,y-1])]
    
    if isempty(indices) == false
        for el in indices

            if el == 1
                if condi1  # Either `v1` is empty or it is already occupied by `ind-1`. A pull move is all but assured.
                    push!(coordinates1,diag1) # Record the possible new coordinates for `ind`.
                    push!(coordinates2,vscoords[1])   # Record the possible new coordinates for `ind-1`.
                elseif condi2 
                    push!(coordinates1,diag1) 
                    push!(coordinates2,vscoords[2])   
                end

            elseif el == 2
                if condi2 
                    push!(coordinates1,diag2) 
                    push!(coordinates2,vscoords[2]) 
                elseif condi3 
                    push!(coordinates1,diag2) 
                    push!(coordinates2,vscoords[3])  
                end

            elseif el == 3 
                if condi3 
                    push!(coordinates1,diag3) 
                    push!(coordinates2,vscoords[3])  
                elseif condi4 
                    push!(coordinates1,diag3) 
                    push!(coordinates2,vscoords[4]) 
                end

            elseif el == 4 
                if condi4 
                    push!(coordinates1,diag4) 
                    push!(coordinates2,vscoords[4])
                elseif condi1 
                    push!(coordinates1,diag4) 
                    push!(coordinates2,vscoords[1])
                end
            end
        end

    end
    


    coordinatesi=[] # Possible coordinates for `ind`
    coordinatesii=[] # Possible coordinates for `ind+1`.
    if isempty(indices2) == false
        for el in indices2

            if el == 1
                if conditio1  # Either `v1` is empty or it is already occupied by `ind+1`. A pull move is all but assured.
                    push!(coordinatesi,diag1) # Record the possible new coordinates for `ind`.
                    push!(coordinatesii,vscoords[1])   # Record the possible new coordinates for `ind+1`.
                elseif conditio2 
                    push!(coordinatesi,diag1) 
                    push!(coordinatesii,vscoords[2])   
                end

            elseif el == 2
                if conditio2 
                    push!(coordinatesi,diag2) 
                    push!(coordinatesii,vscoords[2]) 
                elseif conditio3 
                    push!(coordinatesi,diag2) 
                    push!(coordinatesii,vscoords[3])  
                end

            elseif el == 3 
                if conditio3 
                    push!(coordinatesi,diag3) 
                    push!(coordinatesii,vscoords[3])  
                elseif conditio4 
                    push!(coordinatesi,diag3) 
                    push!(coordinatesii,vscoords[4]) 
                end

            elseif el == 4 
                if conditio4 
                    push!(coordinatesi,diag4) 
                    push!(coordinatesii,vscoords[4])
                elseif conditio1 
                    push!(coordinatesi,diag4) 
                    push!(coordinatesii,vscoords[1])
                end
            end
        end

    end


    # Now I have in `coordinates1,coordinates2` the potential new positions for `ind,ind-1`. The length of these arrays is 
    # the number of possible pull moves por the given index.
    numpull1=length(coordinates1)

    # Now I have in `coordinatesi,coordinatesii` the potential new positions for `ind,ind+1`. The length of these arrays is 
    # the number of possible pull moves por the given index.
    numpull2=length(coordinatesi)
    numpull=numpull1+numpull2

    # I turn the 1 dimensional arrays into matrices for easier use.
    coordinates1=transpose(hcat(coordinates1...))
    coordinates2=transpose(hcat(coordinates2...))

    # I turn the 1 dimensional arrays into matrices for easier use.
    coordinatesi=transpose(hcat(coordinatesi...))
    coordinatesii=transpose(hcat(coordinatesii...))

    return (numpull,coordinates1,coordinates2,coordinatesi,coordinatesii) # Returns everything neccesary to reproduce the final configuration.
end


























#Next, I write a function which counts all of the possible pull moves for a given state.
"""
    countpull2D(N,edo,HPlist)

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; counts all of the possible pull moves. It also outputs the necessary
coordinates to perform one of the listed moves. 
"""
function countpull2D(N,edo,HPlist)

    # Count the moves and store the coordinates for the first monomer.
    t1,coords2,coords1=countFirst2D(N,edo,HPlist)

    # Count the moves and store the coordinates for the last monomer.
    t2,coords4,coords3=countLast2D(N,edo,HPlist)

    


    # Count the moves and store the coordinates for the middle amino acids. 
    t3=0
    tmid=0
    indexb=Int8[] # Stores the value of the middle indices that have non empty arrays of possible backwards pull moves.
    indexf=Int8[] # Stores the value of the middle indices that have non empty arrays of possible forwards pull moves.
    npullpindexb=Int8[] # Stores the number of backwards pull moves for each middle index.
    # Length may be shorter than the number of middle indices given that I only store the coordinates for non empty arrays.
    npullpindexf=Int8[] # Stores the number of forwards pull moves for each middle index.
    middlecoords1=[] # Contains the coordinates for `ind`, where `ind` is an index from the middle of the chain.
    middlecoords2=[] # Contains the coordinates for `ind-1`.
    middlecoords3=[] # Contains the coordinates for `ind`.
    middlecoords4=[] # Contains the coordinates for `ind+1`.
    for j in 2:length(HPlist)-1
        tj,coords5,coords6,coords7,coords8=countMiddle2D(N,j,edo,HPlist)
        # I have all the neccesary information. But first I need to check if the given index
        # posseses at least a possible pull move.
        t3=t3+tj
        if isempty(coords5) == false
            tmid=tmid+size(coords5)[1]
            push!(indexb,j)
            push!(npullpindexb,size(coords5)[1])
            push!(middlecoords1,coords5)
            push!(middlecoords2,coords6)
        end

        if isempty(coords7) == false
            push!(indexf,j)
            push!(npullpindexf,size(coords7)[1])
            push!(middlecoords3,coords7)
            push!(middlecoords4,coords8)
        end
    end

    totalpull=t1+t2+t3

    matrixb=hcat(indexb,npullpindexb)
    matrixf=hcat(indexf,npullpindexf)

    return (totalpull,tmid,matrixb,matrixf,coords1,coords2,coords3,coords4,middlecoords1,middlecoords2,middlecoords3,middlecoords4)
end


















# Before writing the function that performs the actual move, I need a function which returns the index being pullled
# in one of the `middlecoords`matrices from my function `countpull2D`.
"""
    middleInd2D(mp,matrix)

Given a number `mp`, and a matrix `matrx` containing the values of the middle indices with possible pull moves in 
the first column and the number of pull moves for said indices in the second column;returns the corresponding 
index `ind` being moved and the position in one of my `middlecoords[ind][position,:]` matrices.
"""
function middleInd2D(mp,matrix)
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




















# Now that I have the necessary function to count the number of pull moves and the coordinates of such moves, 
# I write a function to actually perform the moves.
"""
    pullMove2D(N,edo,HPlist)

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; chooses an index and performs a full pull move for the generated index. 
"""
function pullMove2D(N,edo,HPlist)
        
    # `newedo` will contain the amino acids´ final positions.
    newedo=copy(edo)

    # Generate the list of possible pull moves.
    totalpull,tmid,matrixb,matrixf,coords1,coords2,coords3,coords4,middlecoords1,middlecoords2,middlecoords3,middlecoords4=countpull2D(N,edo,HPlist)

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
        newedo[ind+1,:]=coords2[m,:] # Record the change.
        indpulled=ind # Function returns the index pulled in the move.
        newcoord1=newedo[ind,:] # Function returns the new coordinate for the pulled index.
        newcoord2=newedo[ind+1,:] 
        dirpulled= forwards


        # Now I have to check whether the current configuration is valid.
        stateconf=validConf(N,ind+1,newedo,HPlist,forwards)
        for k in ind+2:length(HPlist)
            if stateconf == false
                np=edo[k-2,:] # New position.
                newedo[k,:]=np 
                stateconf=validConf(N,k,newedo,HPlist,forwards) # Check whether the new configuration is valid.
            else
                break
            end
        end
        


    elseif  s2 ≤ m ≤ s3
        mp=m-(s2-1) # `mp is the position` for the moves corresponding to the final monomer in the chain.
        ind=length(HPlist)
        newedo[ind,:]=coords3[mp,:] # Record the change.
        newedo[ind-1,:]=coords4[mp,:] # Record the change.
        indpulled=ind # Function returns the index pulled in the move.
        newcoord1=newedo[ind,:] # Function returns the new coordinate for the pulled index.
        newcoord2=newedo[ind-1,:] 
        dirpulled= backwards

        # Now I have to check whether the current configuration is valid.
        stateconf=validConf(N,ind-1,newedo,HPlist,backwards) 
        for k in ind:-1:3
            if stateconf == false
                np=edo[k,:] # New position.
                newedo[k-2,:]=np 
                stateconf=validConf(N,k-2,newedo,HPlist,backwards) # Check whether the new configuration is valid.
            else
                break
            end
        end


    elseif s4 ≤ m ≤ s5
        mp=m-(s4-1)
        indm,ind,pos=middleInd2D(mp,matrixb)
        newedo[indm,:]=middlecoords1[ind][pos,:]
        indpulled=indm # Function returns the index pulled in the move.
        newcoord1=newedo[indm,:] # Function returns the new coordinate for the pulled index. 
        dirpulled= backwards

        if validConf(N,indm,newedo,HPlist,backwards) != true
            newedo[indm-1,:]=middlecoords2[ind][pos,:] # Record the change.
        end
        newcoord2=newedo[indm-1,:]

        # If all went well, a move has been done, now I have to check whether the current configuration is valid.
        stateconf=validConf(N,indm-1,newedo,HPlist,backwards)
        if stateconf == false
            for k in (indm-2):-1:1
                if stateconf == false
                    np=edo[k+2,:] # New position.
                    newedo[k,:]=np 
                    stateconf=validConf(N,k,newedo,HPlist,backwards) # Check whether the new configuration is valid.
                else
                    break
                end
            end
        end
        


    elseif m ≥ s6 
        mp=m-(s6-1)
        indm,ind,pos=middleInd2D(mp,matrixf)
        newedo[indm,:]=middlecoords3[ind][pos,:]
        indpulled=indm # Function returns the index pulled in the move.
        newcoord1=newedo[indm,:] # Function returns the new coordinate for the pulled index. 
        dirpulled= forwards
        

        if validConf(N,indm,newedo,HPlist,forwards) != true
            newedo[indm+1,:]=middlecoords4[ind][pos,:] # Record the change.
        end
        newcoord2=newedo[indm+1,:]


        # If all went well, a move has been performed, now I have to check whether the current configuration is valid.
        stateconf=validConf(N,indm+1,newedo,HPlist,forwards)
        if stateconf ==  false
            for k in (indm+2):length(HPlist)
                if stateconf == false
                    np=edo[k-2,:] # New position.
                    newedo[k,:]=np 
                    stateconf=validConf(N,k,newedo,HPlist,forwards) # Check whether the new configuration is valid.
                else
                    break
                end
            end
        end



    end
    

    return (newedo,totalpull,indpulled,newcoord1,newcoord2,dirpulled) # Returns everything neccesary to reproduce the final configuration and implement
    # the Metropolis algorithm.
end

















