
using StatsBase
using Statistics
using DelimitedFiles




"""
    periodicInd2D(A,indices)

Given a 2D array `A` and a couple of indices, returns the indices for the equivalent array with periodic boundary conditions.
"""
function periodicInd2D(A,indices)
    lx,ly = size(A)
    ix,iy = indices
    Ix = mod1(ix,lx)
    Iy = mod1(iy,ly)
        
    return [Ix,Iy]
end











"""
    periodicArr2D(A,indices)

Given a 2D array `A` and a couple of indices, returns the value of the equivalent array with periodic boundary conditions.
"""
function periodicArr2D(A,indices)
    lx,ly = size(A)
    ix,iy = indices
    Ix = mod1(ix,lx)
    Iy = mod1(iy,ly)
        
    return A[Ix,Iy]
end










# Next, define function to compute the energy of the current state.
"""
    energy(edo)

Given a 2D array `edo`; returns the energy of the configuration.
"""
function energy(edo)

    en = 0 # Stores the value of the energy
    N = size(edo)[1]

    # Iterate over every single place in the lattice. Employ periodic indices.
    for i in 1:N, j in 1:N
        sind = periodicArr2D(edo,[i,j]) # Stores the spin of the lattice site.
        nns = [periodicArr2D(edo,[i-1,j]),periodicArr2D(edo,[i,j+1]),periodicArr2D(edo,[i+1,j]),periodicArr2D(edo,[i,j-1])] # Stores the spins for the 
        # nearest neighbors to the considered spin site.
        localen = sind*(sum(nns)) # Compute the contribution to the energy coming from the spin located at `[i,j]`.
        en = en+localen
    end

    # We have counted each contribution twice. To fix it, divide by 2. Also, since energy is negative, multiply by -1.
    en = (-1/2)*(en)

    return en
end











"""
    nearestneighborSum(edo,inds)

Given a 2D array `edo` and a position `ind`; returns the sum of the nearest neighbors´spins to the given postion.
"""
function nearestneighborSum(edo,inds)

    s = 0 # Stores the value of the sum.
    i,j = inds
    nns = [periodicArr2D(edo,[i-1,j]),periodicArr2D(edo,[i,j+1]),periodicArr2D(edo,[i+1,j]),periodicArr2D(edo,[i,j-1])] # Stores the spins for the nearest neighbors.
    s = sum(nns)
    return s
end












"""
    energyDif(edo,inds)

Given a 2D array `edo` and a position `inds`; returns  the energy difference between the state given 
by `edo`, and a new state where the spin located at `inds` has been flipped .
"""
function energyDif(edo,inds)
    nns = nearestneighborSum(edo,inds) # Nearest neighbor sum.
    difE = Int((2*periodicArr2D(edo,inds))*(nns))
    return difE
end













"""
    histCondtion(localenergies,percent)

Given a dictionary `localenergies` containing the number of times each energy was visited and a percentage `percent`; counts the incidence of each energy. If the
lowest counted energy is `percent` % of the average, return a boolean value `true`.
"""
function histCondtion(localenergies,percent)
    
    ocurrences = values(localenergies)
    hmin = minimum(ocurrences)
    m = mean(ocurrences)
    err = (m-hmin)/m
    val = hmin ≥ m*(percent*1e-2) 

    return (val,err)
end















"""
    ising2DWL(edoin,numlim2)
Given an initial configuration `edo` and a limit number of iterations `numlim2`; returns a dictionary whose keys are the 
visited energies, and it´s values are the logarithms of the energy densities.
"""
function ising2DWL(edo,numlim2)

    N = size(edo)[1]
    mcsweep = N^2 # Monte Carlo sweep/step size.
    sizeSyst = 2*mcsweep
    
    lnf = 1.0 # Initial value of updating value `log(f)`.
    enDensityDict = Dict{Int64,Float64}(-sizeSyst => 0) # Dictionary for the energy densities. Keys are all the possible energy values
    # for the system. Values are the energy density (actually log(g(Eᵢ)) ).

    cont1 = 1
    err = 10.0
    for l in 1:27 #27
        println("cont1= ",cont1)
        cont1 = cont1+1 # Update the counters.
        localenergies = Dict{Int64,Int64}(i => 0 for i in keys(enDensityDict)) # Stores the number of times each energy is visted during the current iteration.
        
        hCond = false
        cont2 = 1
        while (cont2 ≤ numlim2) && (hCond == false)

            for k in 1:(10000*mcsweep) # Perform 10^4 Monte Carlo sweeps.
                
                e1 = Int(energy(edo)) # Initial energy.
                # Generate position to flip.
                iflip = rand(1:N)
                jflip = rand(1:N)
                difE = energyDif(edo,[iflip,jflip]) # Compute the energy difference.
                e2 = Int(e1+difE) # Energy of the new state

                # There´s a chance that the computed energies are not yet in the dictionary. In that case, we add them.
                if e2 ∉ keys(enDensityDict)
                    enDensityDict[e2] = 0
                    localenergies[e2] = 0
                    cont2 = 1
                end
            
                # Next, we need to obtain the density of states for the energies.
                lg1 = enDensityDict[e1]
                lg2 = enDensityDict[e2]
    
                lr = log(rand()) # Generate a random number, then compute it´s log.
                pμν = lg1-lg2 # This quantity detrmines the transition probability.
                if lr ≤ pμν # Accept the flip.
                    edo[iflip,jflip] = (-1)*edo[iflip,jflip] # Flip the spin. 
                    enDensityDict[e2] = lg2+lnf # Update the density of states for the new energy (the log actually).
                    localenergies[e2] = localenergies[e2]+1 # Update the "histogram".
    
                else
                    enDensityDict[e1] = lg1+lnf # Update the density of states for the new energy (the log actually).
                    localenergies[e1] = localenergies[e1]+1  # Update the "histogram".
                end            
    
            end

            hCond, err = histCondtion(localenergies,80)
            if hCond == true 
                println("hCond= ",hCond)
            end

            cont2 = cont2+1
        end
         
        # After the first `while` is completed, we need to change the value of  the modifying constant `f`.
        lnf= (lnf/2) # Update `f`.
    end

     # Next, we need to normalize the dictionary.
     mrest = minimum(values(enDensityDict)) # Choose the minimum density of states.
    
     # Declare a normalized dictionary.
     lngE = Dict{Int64,Float64}(energy => (enDensityDict[energy]-mrest+1) for energy in keys(enDensityDict))
    
    return lngE
end















#=
"""
    ising2DWL_bins(edoin,numlim2)
Given an initial configuration `edo` and a limit number of iterations `numlim2`; returns a dictionary whose keys are the 
visited energies, and it´s values are the logarithms of the energy densities.
"""
function ising2DWL_bins(edo,numlim2)

    N = size(edo)[1]
    mcsweep = N^2 # Monte Carlo sweep/step size.
    sizeSyst = 2*mcsweep
    
    lnf = 1.0 # Initial value of updating value `log(f)`.
    energies = Float64[-sizeSyst] # Energies accesible to the Ising model on the square lattice, of size N^2 and with J=1.
    @show(energies)
    enDensityDict = zeros(Float64,length(energies)) # Dictionary for the energy densities. Values are the energy density (actually log(g(Eᵢ)) ).

    cont1 = 1
    err = 100.0
    for l in 1:27 #27
        println("cont1= ",cont1)
        cont1 = cont1+1 # Update the counters.
        localenergies = zeros(Int64,length(enDensityDict))
        
        hCond = false
        cont2 = 1
        while (cont2 ≤ numlim2) && (hCond == false)

            for k in 1:(10000*mcsweep) # Perform 10^4 Monte Carlo sweeps.
                
                e1 = Int(energy(edo)) # Initial energy.
                
                # Generate position to flip.
                iflip, jflip = rand(1:N,2)
                difE = energyDif(edo,[iflip,jflip]) # Compute the energy difference.
                e2 = Int(e1+difE) # Energy of the new state

                if e2 ∉ energies 
                    push!(energies,e2)
                    push!(enDensityDict,0.0)
                    push!(localenergies,0)
                    cont2 = 1
                end

               # Next, we need to obtain the density of states for the energies.
               ind_e1 = searchsortedlast(energies,e1)
               ind_e2 = searchsortedlast(energies,e2)
               lg1 = enDensityDict[ind_e1]
               lg2 = enDensityDict[ind_e2]
    
                lr = log(rand()) # Generate a random number, then compute it´s log.
                pμν = lg1-lg2 # This quantity detrmines the transition probability.
                if lr ≤ pμν # Accept the flip.
                    edo[iflip,jflip] = (-1)*edo[iflip,jflip] # Flip the spin. 
                    enDensityDict[ind_e2] = lg2+lnf # Update the density of states for the new energy (the log actually).
                    localenergies[ind_e2] = localenergies[ind_e2]+1 # Update the "histogram".
    
                else
                    enDensityDict[ind_e1] = lg1+lnf # Update the density of states for the new energy (the log actually).
                    localenergies[ind_e1] = localenergies[ind_e1]+1  # Update the "histogram".
                end            
    
            end

            hmin = minimum(localenergies)
            m = mean(localenergies)
            err = (m-hmin)/m
            hCond = hmin ≥ m*(0.75) 
            if hCond == true 
                println("hCond= ",hCond)
            end

            cont2 = cont2+1
        end
         
        @show(localenergies)
        @show(err)
        # After the first `while` is completed, we need to change the value of  the modifying constant `f`.
        lnf= (lnf/2) # Update `f`.
    end

    # Next, we need to normalize the dictionary.
    mrest = minimum(enDensityDict) # Choose the minimum density of states.
    
    # Declare a normalized dictionary.
    lngE = Dict{Int64,Float64}(energies[k] => (enDensityDict[k]-mrest+1) for k in 1:length(energies))
    
    return lngE
end
=#











"""
    sExp(E,T,lnge,λ)
Given a temperature `T`, an energy `E`, the logarithm of the enegy density `lnge` and a number `λ`; returns e^{ln(g(e))-βE-λ}.
"""
function sExp(E,T,lnge,λ)
    β = 1/T
    t = (E*β)
    res = exp(lnge-t-λ)
    return res
end
    










"""
    determineλ(lngE,T)
Given a dictionary `lngE` whose keys(values) are energies(log of energy densities), and a temperature; returns the largest value
λ of the differences between energy density `ln(g(E))` and `βE`
"""
function determineλ(lngE,T)
    β = 1/T
    difs = Float64[]
    for el in lngE
        e,lnge = el
        dif = lnge-(β*e)
        push!(difs,dif)
    end
    λ = maximum(difs)
    return λ
end















"""
    ising2DWL_UCFS(edoin,ti,tf,l,numlim2)
 Given an initial state (2D array) `edoin`, an initial (final) temperature `ti` (`tf`), a length for the temperature array `l`
 and a limit number of iterations `numlim2`; returns a matrix containing the thermodynamic variables´ values at each desired temperature.    
"""
function ising2DWL_UCFS(edoin,ti,tf,l,numlim2)

    N = size(edoin)[1] # N^2 is the total number of spins.
    temps = range(ti,stop=tf,length=l) 
    us = ones(Float64,l) # This array will contain the internal energy per spin for each of the temperatures in `temps`.
    cs = ones(Float64,l) 
    Fs = ones(Float64,l)
    Ss = ones(Float64,l)

    lngE = ising2DWL(edoin,numlim2) # Compute the natural logarithm of the energy densities.

    for i in 1:length(temps)
        T = temps[i]
        𝚬 = 0 # Average energy.
        𝚬sq = 0 # Average of the squared energy.
        𝚭 = 0 # Partition function, normalizes the average energy.
        λ = determineλ(lngE,T)
        for element in lngE
            e,lnge = element # Extract the energy and natural logarithm of the energy density.
            z = sExp(e,T,lnge,λ)
            𝚬 = 𝚬+(e*z)
            𝚬sq = 𝚬sq+((e^2)*z)
            𝚭 = 𝚭+z
        end
        
        uT = (𝚬/𝚭) # Internal energy for the given temperature.
        cT = ((𝚬sq/𝚭)-(uT^2))/(T^2) 
        fT = -T*(λ+log(𝚭))
        entropyS = (uT-fT)/T
        
        us[i] = uT/(N^2) # Energy per spin.
        cs[i] = cT/(N^2)
        # Specific heat per spin.
        Fs[i] = fT/(N^2)
        Ss[i] = entropyS/(N^2)
    end

    infoM = ones(Float64,(l,5))
    infoM[:,1] = temps
    infoM[:,2] = us
    infoM[:,3] = cs
    infoM[:,4] = Fs
    infoM[:,5] = Ss

    return (infoM,lngE)
end










# Next, do a little trial run.

#testedo = ones(Float64,(10,10))
#display(testedo)

#res = ising2DWL(testedo,500)


#writedlm("/Users/pedroruiz/Desktop/Diego/PF/isingWL/Data/WLlngE10",res)



