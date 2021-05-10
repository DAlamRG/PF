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
    difes=zeros(Float64,ns) # Stores the energy difference ΔH between states.
    enstates=zeros(Float64,ns+1) # Stores the energy of each state visited during the simulation.
    
    # Compute the energy of the initial state and store it.
    enestados[1]=energy(N,edo,HPlist)




    # Now, apply Metropolis.
    
    #Generate new states by performing pull moves.
    for l in 1:ns 

        # Generate a new state.
        newred,newedo,totalpull = pullMove2D(N,edo,HPlist)
        newenergy=energy(N,newedo,HPlist) # Compute the energy after the pull move.
        
        ΔH=newenergy-enstates[l]
        difes[l]=ΔH
        

        if ΔH ≤ 0 # We accept the configuration change.
            red=copy(newred)
            edo=copy(newedo)
            enstates[l+1]=newenergy # Write the new energy.
        else
            r=rand() # Random number, sampled from uniform distribution over [0,1].
            exponential=exp(-β*ΔH)
            if difE == 4
                if r < as[1] # Se voltea el spin con probabilidad dada por la exponencial
                    red[kx,ky]=(-1)*red[kx,ky]
                    enestados[l+1]=enestados[l]+(difE) 
                    magestados[l+1]=magestados[l]+2*(red[kx,ky]) 
                else
                    enestados[l+1]=enestados[l] # No se voltea el spin y las cantidades permanecen iguales
                    magestados[l+1]=magestados[l]
                end
            elseif difE == 8
                if r < as[2] # Se voltea el spin con probabilidad dada por la exponencial
                    red[kx,ky]=(-1)*red[kx,ky]
                    enestados[l+1]=enestados[l]+(difE) 
                    magestados[l+1]=magestados[l]+2*(red[kx,ky]) 
                else
                    enestados[l+1]=enestados[l] # No se voltea el spin y las cantidades permanecen iguales
                    magestados[l+1]=magestados[l]
                end
            end
        end
        redarray[:,:,l+1]=red[:,:]
    end
    
    return (redarray,enestados)
end


