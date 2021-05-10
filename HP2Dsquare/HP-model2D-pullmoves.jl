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




















@enum directions::Int8 begin
        forwards=1
        backwards=2
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
                elseif dy > 1
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
    dx=x-xplus
    dy=y-yplus
    rel1= dx > 0 && dy == 0
    rel2= dx == 0 && dy < 0
    rel3= dx < 0 && dy == 0
    rel4= dx == 0 && dy > 0 # Conditions single out possible adjacent spaces to move the second amino acid to.

    vs=makecross2D(red,[x,y]) # These are the nearest neighbors to `ind`.
    vscoords=[[x-1,y],[x,y+1],[x+1,y],[x,y-1]]

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
    deleteinds=Int8[]
    for k in 1:length(indices1)
        xa,ya=coordinates2[k] # Possible new coordinates for the first amino acid.
        if validConfEnd2D(red,1,[xa,ya],edo,HPlist) != true
            push!(deleteinds,k) # Store the indices that don´t satisfy the condition.
        end
    end
    deleteat!(coordinates1,deleteinds)
    deleteat!(coordinates2,deleteinds) # Delete possible coordinates.

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
    dx=x-xminus
    dy=y-yminus
    rel1= dx > 0 && dy == 0
    rel2= dx == 0 && dy < 0
    rel3= dx < 0 && dy == 0
    rel4= dx == 0 && dy > 0 # Conditions single out possible adjacent spaces to move the second to last amnino acid to.

    vs=makecross2D(red,[x,y]) # These are the nearest neighbors to `ind`.
    vscoords=[[x-1,y],[x,y+1],[x+1,y],[x,y-1]]

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
            push!(coordinates2,[xv-1,yv])
        elseif el == right
            push!(coordinates2,[xv,yv+1])
        elseif el == down
            push!(coordinates2,[xv+1,yv])
        elseif el == left
            push!(coordinates2,[xv,yv-1])
        end
    end

    # Now I have the coordinates for `ind` and `ind-1`. I proceed to check whether the new configuration would be valid.
    deleteinds=Int8[]
    for k in 1:length(indices1)
        xa,ya=coordinates2[k] # Possible new coordinates for the last amino acid.
        if validConfEnd2D(red,length(HPlist),[xa,ya],edo,HPlist) != true
            push!(deleteinds,k) # STore the indices that don´t satisfy the condition.
        end
    end
    deleteat!(coordinates1,deleteinds)
    deleteat!(coordinates2,deleteinds) # Delete possible coordinates.

    numpull=length(coordinates1)
    # Now we have the possible coordinates for the last two monomers, as well as the number of possible pull moves for the last 
    # amino acid. Fo easier use, I turn the arrays `coordinates1,coordinates2` into a `numpull×2` matrices.

    coordinates1=transpose(hcat(coordinates1...))
    coordinates2=transpose(hcat(coordinates2...))
    return (numpull,coordinates1,coordinates2)
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
    xplus=edo[ind+1,1]
    yplus=edo[ind+1,2] # Position of `ind+1`
    xminus=edo[ind-1,1]
    yminus=edo[ind-1,2] # Position of `ind-1`
    x=edo[ind,1]
    y=edo[ind,2] # Position of `ind`
    
    dx=x-xplus
    dy=y-yplus
    lx=x-xminus
    ly=y-yminus

    rel1= dx > 0 && dy == 0
    rel2= dx == 0 && dy < 0
    rel3= dx < 0 && dy == 0
    rel4= dx == 0 && dy > 0 # Conditions single out possible diagonal spaces to move to.

    rela1= lx > 0 && ly == 0
    rela2= lx == 0 && ly < 0
    rela3= lx < 0 && ly == 0
    rela4= lx == 0 && ly > 0 # Conditions single out possible diagonal spaces to move to.

    # Next we store the possible diagonal spaces.
    diagspaces1=ones(Int8,4)
    if rel1 
        diagspaces1[1]=red[x-1,y+1]
        diagspaces1[4]=red[x-1,y-1]
    elseif rel2 
        diagspaces1[1]=red[x-1,y+1]
        diagspaces1[2]=red[x+1,y+1]
    elseif rel3 
        diagspaces1[2]=red[x+1,y+1]
        diagspaces1[3]=red[x+1,y-1]
    elseif rel4 
        diagspaces1[3]=red[x+1,y-1]
        diagspaces1[4]=red[x-1,y-1]
    end

    diagspaces2=ones(Int8,4)
    if rela1 
        diagspaces2[1]=red[x-1,y+1]
        diagspaces2[4]=red[x-1,y-1]
    elseif rela2 
        diagspaces2[1]=red[x-1,y+1]
        diagspaces2[2]=red[x+1,y+1]
    elseif rela3 
        diagspaces2[2]=red[x+1,y+1]
        diagspaces2[3]=red[x+1,y-1]
    elseif rela4 
        diagspaces2[3]=red[x+1,y-1]
        diagspaces2[4]=red[x-1,y-1]
    end


    v1=red[x-1,y]
    v2=red[x,y+1]
    v3=red[x+1,y]
    v4=red[x,y-1]  # Possible sites to move `ind-1,ind+1` to.

    cond1= edo[ind-1,:] == [x-1,y]
    cond2= edo[ind-1,:] == [x,y+1]
    cond3= edo[ind-1,:] == [x+1,y]
    cond4= edo[ind-1,:] == [x,y-1] # `edo[ind-1,:]` is the current location on the array `red` of `ind-1`.

    condit1= edo[ind+1,:] == [x-1,y]
    condit2= edo[ind+1,:] == [x,y+1]
    condit3= edo[ind+1,:] == [x+1,y]
    condit4= edo[ind+1,:] == [x,y-1] # `edo[ind+1,:]` is the current location on the array `red` of `ind+1`.
    
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
    if isempty(indices) == false
        for el in indices

            if el == 1
                if condi1  # Either `v1` is empty or it is already occupied by `ind-1`. A pull move is all but assured.
                    push!(coordinates1,[x-1,y+1]) # Record the possible new coordinates for `ind`.
                    push!(coordinates2,[x-1,y])   # Record the possible new coordinates for `ind-1`.
                elseif condi2 
                    push!(coordinates1,[x-1,y+1]) 
                    push!(coordinates2,[x,y+1])   
                end

            elseif el == 2
                if condi2 
                    push!(coordinates1,[x+1,y+1]) 
                    push!(coordinates2,[x,y+1]) 
                elseif condi3 
                    push!(coordinates1,[x+1,y+1]) 
                    push!(coordinates2,[x+1,y])  
                end

            elseif el == 3 
                if condi3 
                    push!(coordinates1,[x+1,y-1]) 
                    push!(coordinates2,[x+1,y])  
                elseif condi4 
                    push!(coordinates1,[x+1,y-1]) 
                    push!(coordinates2,[x,y-1]) 
                end

            elseif el == 4 
                if condi4 
                    push!(coordinates1,[x-1,y-1]) 
                    push!(coordinates2,[x,y-1])
                elseif condi1 
                    push!(coordinates1,[x-1,y-1]) 
                    push!(coordinates2,[x-1,y])
                end
            end
        end

    end
    


    coordinatesi=[] # Possible coordinates for `ind`
    coordinatesii=[] # Possible coordinates for `ind+1`.
    if isempty(indices2) == false
        for el in indices2

            if el == 1
                if condi1  # Either `v1` is empty or it is already occupied by `ind+1`. A pull move is all but assured.
                    push!(coordinatesi,[x-1,y+1]) # Record the possible new coordinates for `ind`.
                    push!(coordinatesii,[x-1,y])   # Record the possible new coordinates for `ind+1`.
                elseif condi2 
                    push!(coordinatesi,[x-1,y+1]) 
                    push!(coordinatesii,[x,y+1])   
                end

            elseif el == 2
                if condi2 
                    push!(coordinatesi,[x+1,y+1]) 
                    push!(coordinatesii,[x,y+1]) 
                elseif condi3 
                    push!(coordinatesi,[x+1,y+1]) 
                    push!(coordinatesii,[x+1,y])  
                end

            elseif el == 3 
                if condi3 
                    push!(coordinatesi,[x+1,y-1]) 
                    push!(coordinatesii,[x+1,y])  
                elseif condi4 
                    push!(coordinatesi,[x+1,y-1]) 
                    push!(coordinatesii,[x,y-1]) 
                end

            elseif el == 4 
                if condi4 
                    push!(coordinatesi,[x-1,y-1]) 
                    push!(coordinatesii,[x,y-1])
                elseif condi1 
                    push!(coordinatesi,[x-1,y-1]) 
                    push!(coordinatesii,[x-1,y])
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
    npullpindexb=ones(Int8,length(HPlist)-2) # Stores the number of backwards pull moves for each middle index.
    npullpindexf=ones(Int8,length(HPlist)-2) # Stores the number of forwards pull moves for each middle index.
    middlecoords1=[] # Contains the coordinates for `ind`, where `ind` is an index from the middle of the chain.
    middlecoords2=[] # Contains the coordinates for `ind-1`.
    middlecoords3=[] # Contains the coordinates for `ind`.
    middlecoords4=[] # Contains the coordinates for `ind+1`.
    for j in 2:length(HPlist)-1
        tj,coords5,coords6,coords7,coords8=countMiddle2D(N,j,edo,HPlist)
        push!(middlecoords1,coords5)
        push!(middlecoords2,coords6)
        push!(middlecoords3,coords7)
        push!(middlecoords4,coords8)
        t3=t3+tj
        tmid=tmid+length(coords5[:,1])
        npullpindexb[j-1]=length(coords5[:,1])
        npullpindexf[j-1]=length(coords7[:,1])
    end

    totalpull=t1+t2+t3

    return (totalpull,tmid,npullpindexb,npullpindexf,coords1,coords2,coords3,coords4,middlecoords1,middlecoords2,middlecoords3,middlecoords4)
end


















# Before writing the function that performs the actual move, I need a function which returns the index being pullled
# in one of the `middlecoords`matrices from my function `countpull2D`.
"""
    middleInd2D(mp,npullpindex)

Given a number `mp`, and an array `npullpindex` containing the number of pull moves for each of the middle indices;
returns the corresponding index `ind` being moved and the position in one of my `middlecoords[ind][position,:]` matrices.
"""
function middleInd2D(mp,npullpindex)
    l=length(npullpindex)
    nb=ones(Int8,l)
    nb[1]=npullpindex[1]
    for i in 2:l
        nb[i]=npullpindex[i]+nb[i-1]
    end
    ind=0
    pos=0
    for i in 1:l
        if mp ≤ nb[i]
            ind=i
            break
        end
    end

    if mp ≤ nb[1]
        pos=mp
    else
        for i in 2:l
            if mp ≤ nb[i]
                pos=mp-nb[i-1]
            end
        end
    end
    
    return (ind,pos)
end




















# Now that I have the necessary function to count the number of pull moves and the coordinates of such moves, 
# I write a function to actually perform the moves.
"""
    pullMove2D(N,edo,HPlist)

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of H,P aminoacids `HPlist`; chooses an index and performs a full pull move for the generated index. 
"""
function pullMove2D(N,edo,HPlist)
    
    # I set up the protein within the array.
    red=makeLattice(N,edo,HPlist)
    display(red)
    
    # `newedo` will contain the amino acids´ final positions.
    newedo=copy(edo)

    # Generate the list of possible pull moves.
    totalpull,tmid,npullpindexb,npullpindexf,coords1,coords2,coords3,coords4,middlecoords1,middlecoords2,middlecoords3,middlecoords4=countpull2D(N,edo,HPlist)

    # Randomly choose one of the pull moves.
    m=rand(1:totalpull)

    # MAke subdivisions according to the type of move.
    s1=length(coords1[:,1])
    s2=s1+1
    s3=s1+length(coords3[:,1])
    s4=s3+1
    s5=s3+tmid
    s6=s5+1

    # Type of pull move depends on the value of `m`.
    if m ≤ s1
        ind=1
        singleMove2D(red,edo[ind,:],coords1[m,:]) # Move `ind`.
        newedo[ind,:]=coords1[m,:] # Record the change.
        singleMove2D(red,edo[ind+1,:],coords2[m,:]) # Makes the move.
        newedo[ind+1,:]=coords2[m,:] # Record the change.

        # Now I have to check whether the current configuration is valid.
        stateconf=validConf(N,ind+1,newedo,HPlist,forwards) 
        for k in ind+2:length(HPlist)
            if stateconf == false
                op=edo[k,:] # Old position.
                np=edo[k-2,:] # New position.
                singleMove2D(red,op,np) 
                newedo[k,:]=np 
                stateconf=validConf(N,k,newedo,HPlist,forwards) # Check whether the new configuration is valid.
            else
                break
            end
        end


    elseif  s2 ≤ m ≤ s3
        mp=m-(s2-1)
        ind=length(HPlist)
        singleMove2D(red,edo[ind,:],coords3[mp,:]) # Move `ind`.
        newedo[ind,:]=coords3[mp,:] # Record the change.
        singleMove2D(red,edo[ind-1,:],coords4[mp,:]) # Makes the move.
        newedo[ind-1,:]=coords4[mp,:] # Record the change.

        # Now I have to check whether the current configuration is valid.
        stateconf=validConf(N,ind-1,newedo,HPlist,backwards) 
        for k in ind:-1:3
            if stateconf == false
                op=edo[k-2,:] # Old position.
                np=edo[k,:] # New position.
                singleMove2D(red,op,np) 
                newedo[k-2,:]=np 
                stateconf=validConf(N,k-2,newedo,HPlist,backwards) # Check whether the new configuration is valid.
            else
                break
            end
        end


    elseif s4 ≤ m ≤ s5
        mp=m-(s4-1)
        ind,pos=middleInd2D(mp,npullpindexb)
        singleMove2D(red,edo[ind,:],middlecoords1[ind][pos,:])
        newedo[ind,:]=middlecoords1[ind][pos,:] 
        
        if validConf(N,ind,newedo,HPlist,backwards) != true
            singleMove2D(red,edo[ind-1,:],middlecoords2[ind][pos,:]) # Makes the move.
            newedo[ind-1,:]=middlecoords2[ind][pos,:] # Record the change.
        end

        # If all went well, a move has been done, now I have to check whether the current configuration is valid.
        stateconf=validConf(N,ind-1,newedo,HPlist,backwards) 
        for k in ind:-1:3
            if stateconf == false
                op=edo[k-2,:] # Old position.
                np=edo[k,:] # New position.
                singleMove2D(red,op,np) 
                newedo[k-2,:]=np 
                stateconf=validConf(N,k-2,newedo,HPlist,backwards) # Check whether the new configuration is valid.
            else
                break
            end
        end


    elseif m ≥ s6 
        mp=m-(s6-1)
        ind,pos=middleInd2D(mp,npullpindexf)
        singleMove2D(red,edo[ind,:],middlecoords3[ind][pos,:])
        newedo[ind,:]=middlecoords3[ind][pos,:] 
        
        if validConf(N,ind,newedo,HPlist,forwards) != true
            singleMove2D(red,edo[ind+1,:],middlecoords4[ind][pos,:]) # Makes the move.
            newedo[ind+1,:]=middlecoords4[ind][pos,:] # Record the change.
        end

        # If all went well, a move has been done, now I have to check whether the current configuration is valid.
        stateconf=validConf(N,ind+1,newedo,HPlist,forwards) 
        for k in ind:length(HPlist)-2
            if stateconf == false
                op=edo[k+2,:] # Old position.
                np=edo[k,:] # New position.
                singleMove2D(red,op,np) 
                newedo[k+2,:]=np 
                stateconf=validConf(N,k+2,newedo,HPlist,forwards) # Check whether the new configuration is valid.
            else
                break
            end
        end


    end

    return (red,newedo,totalpull) # Returns everything neccesary to reproduce the final configuration.
end

















