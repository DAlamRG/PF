# First attempt at programming the Metropolis algorithm for the HP model.

include("./HP-model2D-pullmoves.jl")




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
























# Next, I define  a function which takes the configuration of the protein within a 2Darray and outputs its energy.
"""
    energy(N,edo,HPlist)

Given a 2D array size `N`, a matrix encoding the aminoacids positions `edo`, and an array containing the sequence of H,P aminoacids `HPlist`; 
returns the energy of the configuration.
"""
function energy(N,edo,HPlist)
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
    en=((enp)-2*countH(HPlist))/2 # The total energy is obtained by substracting the number of covalent bonds to `enp` 
    #and dividing the resulting number by two.
    return en
end
























# Now that I am able to generate a new state by performing a pull move, and knowing how to compute the energy of any
# given configuration, I proceed to implement a simulation using the Metropolis-Hastings algorithm.

"""
    HP2Dmet(N,ns,T,edo,HPlist) 

Given a 2D array size `N`, a number of steps `ns`, the temperature T, a matrix encoding the aminoacids positions `edo`, 
and an array containing the sequence of H,P aminoacids `HPlist`; returns the final configuration after performing a simulation 
using the Metropolis-Hastings algorithm.
"""
function HP2Dmet(N,ns,T,edo,HPlist) 
    

    β=1/T
    red=makeLattice(N,edo,HPlist)
    
    redarray=ones(Int8,(N,N,ns+1)) # Multidimensional array whose entries `redarray[:,:,j]` are the matrices`red` at a time t=j.
    redarray[:,:,1]=red # The first array is filled by the input state.

    states=ones(Int8,(length(HPlist),2,ns+1)) # `states` stores the amino acids´ coordinates for each visited configuration.
    states[:,:,1]=edo # First coordinates are those of the input state.
    
    difes=zeros(Float64,ns) # Stores the energy difference ΔH between states.
    enstates=zeros(Float64,ns+1) # Stores the energy of each state visited during the simulation.
    enstates[1]=energy(N,edo,HPlist) # Compute the energy of the initial state and store it.
    
    npullstates=zeros(Float64,ns+1) # Stores the number of pull moves for each configuration.
    npullstates[1]=countpull2D(N,edo,HPlist)[1] # Count the number of pull moves for the initial state and store it.

    
    

    # Now, apply Metropolis-Hastings.
    
    #Generate new states by performing pull moves.
    for l in 2:ns+1 

        # Generate a new state.
        newred,newedo,totalpull = pullMove2D(N,states[:,:,l-1],HPlist)
        newenergy=energy(N,newedo,HPlist) # Compute the energy after the pull move.
        
        ΔH=newenergy-enstates[l-1] # Compute the energy difference.
        
        

        if ΔH ≤ 0 # We accept the configuration change.
            redarray[:,:,l]=newred
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
                redarray[:,:,l]=newred
                states[:,:,l]=newedo
                enstates[l]=newenergy 
                npullstates[l]=totalpull
                difes[l-1]=ΔH
            else # MOve is not accepted, the protein´s configuration stays the same.
                redarray[:,:,l]=redarray[:,:,l-1]
                states[:,:,l]=states[:,:,l-1]
                enstates[l]=enstates[l-1]
                npullstates[l]=npullstates[l-1]
                difes[l-1]=0
            end
        end
    end


    # After the above for loop has ended, I have the final configuration for the protein. I pass the configuration
    # as the function´s output.
    
    return (redarray,states,enstates)
end


#[[20 10];[20 11];[20 12];[20 13];[20 14];[20 15];[20 16];[20 17];[20 18];[20 19];[20 20];[20 21];[20 22];[20 23];[20 24];[20 25];[20 26];[20 27];[20 28];[20 29]]
#[-1,1,-1,1,1,-1,-1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,-1]

info=HP2Dmet(40,40000,0.2,[[20 10];[20 11];[20 12];[20 13];[20 14];[20 15];[20 16];[20 17];[20 18];[20 19];[20 20];[20 21];[20 22];[20 23];[20 24];[20 25];[20 26];[20 27];[20 28];[20 29]],[-1,1,-1,1,1,-1,-1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,-1]