
using StatsBase
using Statistics
# I try my hand at the Ising model using the Wang-Landau algorithm.



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










# Next, define function to compute the energy of the current state.
"""
    energy(edo)

Given a 2D array `edo`; returns the energy of the configuration.
"""
function energy(edo)

    en=0 # Stores the value of the energy
    N=size(edo)[1]

    # Iterate over every single place in the lattice. Employ periodic indices.
    for i in 1:N, j in 1:N
        sind=periodicArr2D(edo,[i,j]) # Stores the spin of the lattice site.
        nns=[periodicArr2D(edo,[i-1,j]),periodicArr2D(edo,[i,j+1]),periodicArr2D(edo,[i+1,j]),periodicArr2D(edo,[i,-1j])] # Stores the spins for the 
        # nearest neighbors to the considered spin site.
        localen=sind*(sum(nns)) # COmpute the contribution to the energy coming from the spin located at `[i,j]`.
        en=en+localen
    end

    # We have counted each contribution twice. To fix it, divide by 2. ALso, since energy is negative, multiply by -1.
    en=(-1/2)*(en)

    return en
end











"""
    nearestneighborSum(edo,inds)

Given a 2D array `edo` and a position `ind`; returns the sum of the nearest neighbors´spins to the given postion.
"""
function nearestneighborSum(edo,inds)

    s=0 # Stores the value of the energy
    i,j=inds
    nns=[periodicArr2D(edo,[i-1,j]),periodicArr2D(edo,[i,j+1]),periodicArr2D(edo,[i+1,j]),periodicArr2D(edo,[i,-1j])] # Stores the spins for the nearest neighbors.
    s=sum(nns)
    return s
end












"""
    energyDif(edo,inds)

Given a 2D array `edo` and a position `inds`; returns  the energy difference between the state given 
by `edo`, and a new state where the spin located at `inds` has been flipped .
"""
function energyDif(edo,inds)
    nns=nearestneighborSum(edo,inds) # Nearest neighbor sum.
    difE=Int((2*periodicArr2D(edo,inds))*(nns))
    return difE
end













"""
    histCondtion(localenergies)

Given a dictionary `localenergies` containing the number of times each energy was visited; counts the incidence of each energy. If the
lowest counted energy is 80 % of the average, return a boolean value `true`.
"""
function histCondtion(localenergies)
    val = true
    ocurrences=values(localenergies)
    m1=minimum(ocurrences)
    μ=mean(ocurrences)
    μcomp=(0.5)*μ

    if m1 ≤ μcomp
        val = false
    end
    return val
end












"""
    ising2DWL(edoin)
"""
function ising2DWL(edo)

    N=size(edo)[1]
    mcsweep=N^2 # Monte Carlo sweep/step size.
    numlim1=10^(4) # Limits the number of times the histogram is updated. Avoids non-ending algorithm.
    numlim2=10^(5)
    
    lnf=1 # Initial value of updating value `f`.
    enlim=Int(2*(N^2)) # Upper energy limit for the N^2 spin Ising system.
    enDensityDict=Dict(i => log(exp(1)) for i in -enlim:enlim) # Dictionary for the energy densities. Keys are all the possible energy values
    # for the system. Values are the energy density (actually log(g(Eᵢ)) ).



    cont1=1
    while (cont1 ≤ numlim1) && (lnf ≥ 10^(-4)) # Stop when modification factor is small enough.
        cont1=cont1+1 # Update the counter.
        println("cont1= ",cont1)
        localenergies=Dict(i => 0 for i in -enlim:enlim) # Stores the number of times each energy is visted during the current iteration.
        
        hCond=false
        cont2=1
        while (cont2 ≤ numlim2) && (hCond == false)
            cont2=cont2+1
            

            for k in 1:(100*mcsweep) # Perform 10^3 Monte Carlo sweeps.
                
                e1=Int(energy(edo)) # Initial energy.
                # Generate position to flip.
                iflip=rand(1:N)
                jflip=rand(1:N)
                difE=energyDif(edo,[iflip,jflip]) # Compute the energy difference.
                e2=Int(e1+difE) # Energy of the new state
    
                # Next, we need to obtain the density of states for the energies.
                lg1=enDensityDict[e1]
                lg2=enDensityDict[e2]
    
                r=rand(0:1) # Generate a random number.
                if r ≤ exp(lg1-lg2) # Accept the flip.
                    edo[iflip,jflip]=-edo[iflip,jflip] # Flip the spin. 
                    enDensityDict[e2]=lg2+lnf # Update the density of states for the new energy (the log actually).
                    localenergies[e2]=localenergies[e2]+1 # Update the "histogram".
    
                else
                    enDensityDict[e1]=lg1+lnf # Update the density of states for the new energy (the log actually).
                    localenergies[e1]=localenergies[e1]+1  # Update the "histogram".
                end            
    
            end

            hCond=histCondtion(localenergies)
        end

        display(localenergies)
            
        # After the first `while` is completed, we need to change the value of  the modifying constant `f`.
        lnf= (1/2)*lnf # Update `f`.
        localenergies=Dict(i => 0 for i in -enlim:enlim) # Reset histogram.
    end

    return enDensityDict
end










"""
    sExp(E,T)
Given a temperature `T` and an energy `E`; returns e^{-βE}.
"""
function sExp(E,T)
    β=1/T
    t=-(E*β)
    return exp(t)
end
    









"""
    ising2DWL_U(edoin,ti,tf,l)
 Given an initial state (2D array) `edoin`, an initial (final) temperature `ti` (`tf`) and a length for the temperature array `l`; 
 returns an array containing the internal energy per spin for each of the temperatures in the range.    
"""
function ising2DWL_U(edoin,ti,tf,l)

    N=size(edoin)[1]
    temps=range(ti,stop=tf,length=l)
    us=ones(Float64,l) # This array will contain the internal energy per spin for each of the temperatures in `temps`.

    lngE=ising2DWL(edoin) # Compute the natural logarithm of the energy densities.

    for i in 1:length(temps)
        T=temps[i]
        s1=0
        s2=0
        for el in lngE
            E,lnge=el # Extract the energy and natural logarithm of the energy density.
            ge=exp(lnge) # Compute the energy density corresponding to the enegy `E`.
            common=sExp(E,T)*ge
            s1=s1+(E*common)
            s2=s2+common
        end
        uT=(s1/s2) # INternal energy for the given temperature.
        us[i]=uT/(N^2) # Energy per spin.
    end

    return us
end






# Next, do a little trial run.

testedo=Int64[[1 -1 ];[-1 1]]


