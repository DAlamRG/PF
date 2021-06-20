# First attempt at programming the Metropolis algorithm for the HP model.

include("./HP-model2D-pullmoves-triangle.jl")




# First, I need a function that, given a sequence of aminoacids, returns the number of covalent bonds.
"""
    countH(HPlist)

Given a vector (H,H,P,H,P,...)=(-1,-1,1,-1,1,...), counts the number of H-H bonds within the array.
"""
function countH(HPlist)
    s=0
    for i in 1:(length(HPlist)-1)
        if HPlist[i] == 1
            continue
        elseif HPlist[i] == -1
            if HPlist[i+1] == 1
                continue
            elseif HPlist[i+1] == -1
                s=s+1
            end
        end
    end
    return s
end
























# Next, I define  a function which takes the configuration of the protein within a 2Darray and outputs its energy. It also takes the geometry of
# the lattice into account.
"""
    energy(N,edo,HPlist,geometry)

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`, an array containing the sequence of H,P aminoacids `HPlist`,
and the geometry of the lattice; returns the energy of the configuration.
"""
function energy(N,edo,HPlist,geometry)
    red=makeLattice(N,edo,HPlist)
    
    # I create an auxiliary array with periodic boundary conditions.
    redaux=zeros(Int8,(N+2,N+2))
    redaux[2:end-1,2:end-1]=red
    redaux[end,2:end-1]=red[1,:]
    redaux[1,2:end-1]=red[end,:]
    redaux[2:end-1,end]=red[:,1]
    redaux[2:end-1,1]=red[:,end]
    
    # Next, I iterate over the protein´s vertices and compute the number of H-H bonds.
    # To avoid counting the covalent bonds, I substract two times the value of `countBH(HPlist)`.
    enp=0

    if geometry == square2D
        for i in 1:length(HPlist)
            if HPlist[i]==1 
                continue
            elseif HPlist[i]==-1
                x=edo[i,1]+1
                y=edo[i,2]+1 # I have to shift both indices by one place because we are working on the extended array `redaux`.
                cruz=[redaux[x+1,y],redaux[x,y+1],redaux[x-1,y],redaux[x,y-1]] # I define an array with the next closest neighbors.
                c=count(i->(i==-1),cruz)
                enp=enp+(c)
            end
        end

    elseif geometry == triangular2D
        for i in 1:length(HPlist)
            if HPlist[i]==1 
                continue
            elseif HPlist[i]==-1
                x=edo[i,1]+1
                y=edo[i,2]+1 # I have to shift both indices by one place because we are working on the extended array `redaux`.
                # I define an array with the nearest neighbors.
                nn=[redaux[x-1,y],redaux[x,y+1],redaux[x+1,y+1],redaux[x+1,y],redaux[x,y-1],redaux[x-1,y-1]] 
                c=count(i->(i==-1),nn)
                enp=enp+(c)
            end
        end
    end



    en=((enp)-2*countH(HPlist))/2 # The total energy is obtained by substracting the number of covalent bonds to `enp` 
    #and dividing the resulting number by two.
    return -en
end

























# Now that I am able to generate a new state by performing a pull move, and knowing how to compute the energy of any
# given configuration, I proceed to implement a simulation using the Metropolis-Hastings algorithm.

"""
    HP2Dmet(N,nums,T,protein) 

Given a 2D array size `N`, a number of Monte-Carlo sweeps `nums`, the temperature T, and a structure `protein` encoding the protein´s configuration; 
returns the final configuration and the visited energies after performing a simulation using the Metropolis-Hastings algorithm.
"""
function HP2Dmet(N,nums,T,protein) 
    

    β=1/T
    edo=protein.edo
    HPlist=protein.HPlist
    geometry=protein.geometry
    ns=nums*length(HPlist) # Total number of iterations

    
    states=ones(Int8,(length(HPlist),2,ns+1)) # `states` stores the amino acids´ coordinates for each visited configuration.
    states[:,:,1]=edo # First coordinates are those of the input state.
    
    difes=zeros(Float64,ns) # Stores the energy difference ΔH between states.
    enstates=zeros(Float64,ns+1) # Stores the energy of each state visited during the simulation.
    enstates[1]=energy(N,edo,HPlist,geometry) # Compute the energy of the initial state and store it.
    
    npullstates=zeros(Float64,ns+1) # Stores the number of pull moves for each configuration.
    npullstates[1]=totalpull2Dg(N,edo,HPlist,geometry)[1] # Count the number of pull moves for the initial state and store it.


    

    # Now, apply Metropolis-Hastings. The function´s output dependds on the geometry.
    if geometry == square2D
        pulledindices=zeros(Int8,ns) # Stores the moved indices. If no index is moved, the value is zero.
        newcoords1=zeros(Int8,(ns,2)) # Stores the new coordinate for the pulled index. 
        newcoords2=zeros(Int8,(ns,2)) # Stores the new coordinate for the momomer directly before/after the one being moved. 
        dirs=directions[] # Store the direction in which the chain was pulled. If no pull move was performed, the direction is `nonetaken`.

        #Generate new states by performing pull moves.
        for l in 2:ns+1 
            # Generate a new state.
            newedo,totalpull,indpulled,newcoord1,newcoord2,dirpulled = pullMove2DGeneral(N,states[:,:,l-1],HPlist,geometry)
            newenergy=energy(N,newedo,HPlist,geometry) # Compute the energy after the pull move.
            
            ΔH=newenergy-enstates[l-1] # Compute the energy difference.
            
    
            if ΔH ≤ 0 # We accept the configuration change.
                states[:,:,l]=newedo
                enstates[l]=newenergy 
                npullstates[l]=totalpull
                difes[l-1]=ΔH
                pulledindices[l-1]=indpulled
                newcoords1[l-1,:]=newcoord1
                newcoords2[l-1,:]=newcoord2
                push!(dirs,dirpulled)
    
    
            else
                r=rand() # Random number, sampled from uniform distribution over [0,1].
                exponential=exp(-β*ΔH)
                q=(npullstates[l-1])/(totalpull) # This is the quotient of number of pull moves.
                pστ=q*exponential # This determines the probability of the move in the case ΔH > 0.
    
                if r < pστ # Accept the move.
                    states[:,:,l]=newedo
                    enstates[l]=newenergy 
                    npullstates[l]=totalpull
                    difes[l-1]=ΔH
                    pulledindices[l-1]=indpulled
                    newcoords1[l-1,:]=newcoord1
                    newcoords2[l-1,:]=newcoord2
                    push!(dirs,dirpulled)
                else # Move is not accepted, the protein´s configuration stays the same.
                    states[:,:,l]=copy(states[:,:,l-1])
                    enstates[l]=copy(enstates[l-1])
                    npullstates[l]=copy(npullstates[l-1])
                    difes[l-1]=0
                    pulledindices[l-1]=0
                    newcoords1[l-1,:]=[0,0]
                    newcoords2[l-1,:]=[0,0]
                    push!(dirs,nonetaken)
                end
            end
        end

        newcoords=(newcoords1,newcoords2)
        return (pulledindices,dirs,newcoords)


    elseif geometry == triangular2D
        pulledindices=zeros(Int8,ns) # Stores the moved indices. If no index is moved, the value is zero.
        newcoords=zeros(Int8,(ns,2)) # Stores the new coordinate for the pulled index. 
        dirs=directions[] # Store the direction in which the chain was pulled. If no pull move was peroformed, the direction is `nonetaken`.

        #Generate new states by performing pull moves.
        for l in 2:ns+1 
            
            # Generate a new state.
            newedo,totalpull,indpulled,newcoord,dirpulled = pullMove2DGeneral(N,states[:,:,l-1],HPlist,geometry)
            newenergy=energy(N,newedo,HPlist,geometry) # Compute the energy after the pull move.
            
            ΔH=newenergy-enstates[l-1] # Compute the energy difference.
            
    
            if ΔH ≤ 0 # We accept the configuration change.
                states[:,:,l]=newedo
                enstates[l]=newenergy 
                npullstates[l]=totalpull
                difes[l-1]=ΔH
                pulledindices[l-1]=indpulled
                newcoords[l-1,:]=newcoord
                push!(dirs,dirpulled)

    
            else
                r=rand() # Random number, sampled from uniform distribution over [0,1].
                exponential=exp(-β*ΔH)
                q=(npullstates[l-1])/(totalpull) # This is the quotient of number of pull moves.
                pστ=q*exponential # This determines the probability of the move in the case ΔH > 0.
    
                if r < pστ # Accept the move.
                    states[:,:,l]=newedo
                    enstates[l]=newenergy 
                    npullstates[l]=totalpull
                    difes[l-1]=ΔH
                    pulledindices[l-1]=indpulled
                    newcoords[l-1,:]=newcoord
                    push!(dirs,dirpulled)
                else # Move is not accepted, the protein´s configuration stays the same.
                    states[:,:,l]=copy(states[:,:,l-1])
                    enstates[l]=copy(enstates[l-1])
                    npullstates[l]=copy(npullstates[l-1])
                    difes[l-1]=0
                    pulledindices[l-1]=0
                    newcoords[l-1,:]=[0,0]
                    push!(dirs,nonetaken)
                end
            end
        end

        return (pulledindices,dirs,newcoords)
    end
end



 












# Next, I write a function which performs multiple simulations over an array of temperatures.
"""
    mainHP2Dmet(N,nums,ti,tf,nTs,protein)
Given a 2D array size `N`, a number of Monte-Carlo sweeps `nums` per temperature, an initial(final) temperature `ti(tf)`, the number of 
temperatures to be visited `nTs`, and a structure `protein` encoding the protein´s configuration; returns the final configuration and 
the visited energies after performing a simulation using the Metropolis-Hastings algorithm.
"""
function mainHP2Dmet(N,nums,ti,tf,nTs,protein)

    temperatures=range(ti,stop=tf,length=nTs) # Decalre a range of temperatures to be visited.
    geom=protein.geometry

    if geom == square2D
        
        # Next, I need three arrays which will store all of the `pulledindices,dirs,newcoords`for all of the temperatures. 
        pulledindicesT=Vector{Int8}[] 
        newcoordsT=Tuple[] 
        dirsT=Vector{directions}[]

        # Fill the first entry of the above arrays.
        datatemp1=HP2Dmet(N,nums,temperatures[end],protein)
        push!(pulledindicesT,datatemp1[1])
        push!(dirsT,datatemp1[2])
        push!(newcoordsT,datatemp1[3])

        laststate=reconstructStates2D(N,protein.edo,protein.HPlist,pulledindicesT[1],dirsT[1],newcoordsT[1],protein.geometry)[:,:,end]

        # For each of the remaining temperatures, employ the Metropolis-Hastings algorithm to store the information about the
        # visited configurations at the current temperature.
        for k in 2:length(temperatures)
            temp=reverse(temperatures)[end-(k-1)] # temperature at which the simulation is performed
            proteintemp=Protein2D(laststate,protein.HPlist,protein.geometry) # Protein structure for the simulation.
            pulledindices,dirs,newcoords=HP2Dmet(N,nums,temp,proteintemp)  # Perfomr the simulation.
            push!(pulledindicesT,pulledindices) # Record the results.
            push!(dirsT,dirs)
            push!(newcoordsT,newcoords)
            laststate=reconstructStates2D(N,laststate,protein.HPlist,pulledindicesT[k],dirsT[k],newcoordsT[k],protein.geometry)[:,:,end] # Reconstruct the last configuration from the previous temperature.
        end

        return (pulledindicesT,newcoordsT,dirsT,laststate)







    elseif geom == triangular2D

        # Next, I need three arrays which will store all of the `pulledindices,dirs,newcoords`for all of the temperatures. 
        pulledindicesT=Vector{Int8}[] 
        newcoordsT=Matrix{Int64}[] 
        dirsT=Vector{directions}[]

        # Fill the first entry of the above arrays.
        datatemp1=HP2Dmet(N,nums,temperatures[end],protein)
        push!(pulledindicesT,datatemp1[1])
        push!(dirsT,datatemp1[2])
        push!(newcoordsT,datatemp1[3])

        laststate=reconstructStates2D(N,protein.edo,protein.HPlist,pulledindicesT[1],dirsT[1],newcoordsT[1],protein.geometry)[:,:,end]

        # For each of the remaining temperatures, employ the Metropolis-Hastings algorithm to store the information about the
        # visited configurations at the current temperature.
        for k in 2:length(temperatures)
            temp=reverse(temperatures)[end-(k-1)] # temperature at which the simulation is performed
            proteintemp=Protein2D(laststate,protein.HPlist,protein.geometry) # Protein structure for the simulation.
            pulledindices,dirs,newcoords=HP2Dmet(N,nums,temp,proteintemp)  # Perfomr the simulation.
            push!(pulledindicesT,pulledindices) # Record the results.
            push!(dirsT,dirs)
            push!(newcoordsT,newcoords)
            laststate=reconstructStates2D(N,laststate,protein.HPlist,pulledindicesT[k],dirsT[k],newcoordsT[k],protein.geometry)[:,:,end] # Reconstruct the last configuration from the previous temperature.
        end

        return (pulledindicesT,newcoordsT,dirsT,laststate)
    
    end

end