# First attempt at programming the Metropolis algorithm for the HP model.

include("/Users/pedroruiz/Desktop/Diego/PF/HP3D/HP-model3D-pullmoves.jl")
include("/Users/pedroruiz/Desktop/Diego/PF/HP3D/visualizeHP3D.jl")




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

Given a 3D array size `N`, a matrix encoding the aminoacids positions `edo`, an array containing the sequence of H,P aminoacids `HPlist`,
and the geometry of the lattice; returns the energy of the configuration.
"""
function energy(N,edo,HPlist,geometry)
    red=makeLattice3D(N,edo,HPlist)
    
  
    # Next, I iterate over the protein´s vertices and compute the number of H-H bonds.
    # To avoid counting the covalent bonds, I substract two times the value of `countBH(HPlist)`.
    enp=0

    if geometry == cubic
        for i in 1:length(HPlist)
            if HPlist[i]==1 
                continue
            elseif HPlist[i]==-1
                # I define an array with the nearest neighbors.
                nn=nearestNeighbors3D(red,edo[i,:],geometry)
                c=count(i->(i==-1),nn)
                enp=enp+(c)
            end
        end

    elseif geometry == hcp
        for i in 1:length(HPlist)
            if HPlist[i]==1 
                continue
            elseif HPlist[i]==-1
                # I define an array with the nearest neighbors.
                nn=nearestNeighbors3D(red,edo[i,:],geometry)
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
    HP3Dmet(N,nums,T,protein) 

Given a 3D array size `N`, a number of Monte-Carlo sweeps `nums`, the temperature T, and a structure `protein` encoding the protein´s configuration; 
returns the final configuration and the visited energies after performing a simulation using the Metropolis-Hastings algorithm.
"""
function HP3Dmet(N,nums,T,protein) 
    

    β=1/T
    edo=protein.edo
    HPlist=protein.HPlist
    geometry=protein.geometry
    ns=nums*length(HPlist) # Total number of iterations

    red=makeLattice3D(N,edo,HPlist)
    
    redarray=ones(Int8,(N,N,N,ns+1)) # Multidimensional array whose entries `redarray[:,:,j]` are the matrices`red` at a time t=j.
    redarray[:,:,:,1]=red # The first array is filled by the input state.

    states=ones(Int8,(length(HPlist),3,ns+1)) # `states` stores the amino acids´ coordinates for each visited configuration.
    states[:,:,1]=edo # First coordinates are those of the input state.
    
    difes=zeros(Float64,ns) # Stores the energy difference ΔH between states.
    enstates=zeros(Float64,ns+1) # Stores the energy of each state visited during the simulation.
    enstates[1]=energy(N,edo,HPlist,geometry) # Compute the energy of the initial state and store it.
    
    npullstates=zeros(Float64,ns+1) # Stores the number of pull moves for each configuration.
    npullstates[1]=countpull3D(N,edo,HPlist,geometry)[1] # Count the number of pull moves for the initial state and store it.

    
    

    # Now, apply Metropolis-Hastings.
    
    #Generate new states by performing pull moves.
    for l in 2:ns+1 
        

        # Generate a new state.
        newred,newedo,totalpull = pullMove3D(N,states[:,:,l-1],HPlist,geometry)
        newenergy=energy(N,newedo,HPlist,geometry) # Compute the energy after the pull move.
        
        ΔH=newenergy-enstates[l-1] # Compute the energy difference.
        
        

        if ΔH ≤ 0 # We accept the configuration change.
            redarray[:,:,:,l]=newred
            states[:,:,l]=newedo
            enstates[l]=newenergy 
            npullstates[l]=totalpull
            difes[l-1]=ΔH


        else
            r=rand() # Random number, sampled from uniform distribution over [0,1].
            exponential=exp(-β*ΔH)
            q=(npullstates[l-1])/(totalpull) # This is the quotient of number of pull moves.
            pστ=q*exponential # This determines the probability of the move in the case ΔH > 0.

            if r < pστ # Accept the move.
                redarray[:,:,:,l]=newred
                states[:,:,l]=newedo
                enstates[l]=newenergy 
                npullstates[l]=totalpull
                difes[l-1]=ΔH
            else # Move is not accepted, the protein´s configuration stays the same.
                redarray[:,:,:,l]=copy(redarray[:,:,:,l-1])
                states[:,:,l]=copy(states[:,:,l-1])
                enstates[l]=copy(enstates[l-1])
                npullstates[l]=copy(npullstates[l-1])
                difes[l-1]=0
            end
        end
    end


    # After the above for loop has ended, I have the final configuration for the protein. I pass the configuration
    # as the function´s output.
    
    return (redarray,states,enstates)
end



 




#testproteinT=Protein([[7 2 9];[7 3 9];[7 4 9];[7 5 9];[7 6 9];[7 7 9];[7 8 9];[7 9 9];[7 10 9];[7 11 9]],[1,-1,-1,1,-1,1,-1,1,-1,-1],hcp)

# infoT=HP3Dmet(20,40000,0.3,testproteinT)





# [[7 3];[7 4];[7 5];[7 6];[7 7];[6 7];[6 8];[7 9];[7 10];[7 11]]
