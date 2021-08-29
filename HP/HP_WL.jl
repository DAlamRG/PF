using Base: Float64, String
# This script contains the code to perform a Wang-Landau simulation for the HP and extended models of
# protein folding.

using StatsBase
using Statistics
using DelimitedFiles

include("./Energy.jl")







"""
    energyDif(N,edoi,edof,HPlist,geometry,pfmodel::PF_model)

Given lattice size `N`, an initial/final configuration `edoi/edof`, a list containing the protein's aminoacids `HPlist`,
a geometry `geometry` and an interaction model `pfmodel`; returns the energy difference between the intial and final configurations.
"""
function energyDif(N,edoi,edof,HPlist,geometry,pfmodel::PF_model)
    ΔE = energy(N,edoi,HPlist,geometry,pfmodel)-energy(N,edof,HPlist,geometry,pfmodel)
    return ΔE
end











"""
    histCondtion(localenergies)

Given a dictionary `localenergies` containing the number of times each energy was visited; counts the incidence of each energy. If the
lowest counted energy is 80 % of the average, return a boolean value `true`.
"""
function histCondtion(localenergies)
    val = false
    
    ocurrences = values(localenergies)
    hmin = minimum(ocurrences)
    hmax = maximum(ocurrences)
    c = (hmax-hmin)/(hmax+hmin)

    if c < 0.2
        val = true
    end
    return val
end












"""
    wang_landau(N,edo,HPlist,geometry,numlim2,pfmodel::PF_model)
Given a lattice size `N`, an initial configuration `edo`, an aminoacid list `HPlist`, a geometry `geometry`, a limit number of 
iterations `numlim2` and an interaction model `pfmodel`; returns a dictionary whose keys are the visited energies, and it´s 
values are the logarithms of the energy densities.
"""
function wang_landau(N,edo,HPlist,geometry,numlim2,pfmodel::PF_model)

    mcsweep = length(edo) # Monte Carlo sweep size.
    
    lnf = 1 # Initial value of updating value `log(f)`.
    enDensityDict = Dict{Float64,Float64}(0 => 0) # Dictionary for the energy densities. Keys are all the possible energy values
    # for the system (I start with just one since the range of possible energies depends entirely on the aminoacid list/geometry/interaction model).
    # Values are the energy density (actually log(g(Eᵢ)) ), where g(Eᵢ) is the energy density.

    cont1 = 1
    for l in 1:27 # This is the number of iterations it takes to make lnf sufficiently small.
        cont1 = cont1+1 # Update the counter.
        println("cont1= ",cont1)
        localenergies=Dict{Float64,Int64}() # Stores the number of times each energy is visted during the current iteration. This is my "histogram".
        
        hCond = false # True if the histogram is flat enough.
        cont2 = 1
        while (cont2 ≤ numlim2) && (hCond == false)

            for k in 1:(30*mcsweep) # Perform 10^3 Monte Carlo sweeps.
                
                e1 = energy(N,edo,HPlist,geometry,pfmodel) # Initial energy.
                # Generate new state by performing a pull-move.
                edoaux = pullMove(N,edo,HPlist,geometry)[1]
                e2 = energy(N,edoaux,HPlist,geometry,pfmodel) # Energy of the new state

                # There´s a chance that the computed energies are not yet in the dictionary. In that case, we add them.
                if e1 ∉ keys(enDensityDict)
                    enDensityDict[e1] = 0
                end
                if e2 ∉ keys(enDensityDict)
                    enDensityDict[e2] = 0
                end

                if e1 ∉ keys(localenergies)
                    localenergies[e1] = 0
                end
                if e2 ∉ keys(localenergies)
                    localenergies[e2] = 0
                end
    
                # Next, we need to obtain the density of states for the energies.
                lg1 = enDensityDict[e1]
                lg2 = enDensityDict[e2]
    
                r = rand() # Generate a random number.
                pμν = exp(lg1-lg2)
                if r ≤ pμν # Accept the flip.
                    edo = edoaux # Accept the new configuration.  
                    enDensityDict[e2] = lg2+lnf # Update the density of states for the new energy (the log actually).
                    localenergies[e2] = localenergies[e2] + 1 # Update the "histogram".
    
                else
                    enDensityDict[e1] = lg1+lnf # Update the density of states for the new energy (the log actually).
                    localenergies[e1] = localenergies[e1] + 1  # Update the "histogram".
                end            
    
            end

            hCond = histCondtion(localenergies)
            if hCond == true 
                println("hCond = ",hCond)
            end

            cont2 = cont2 + 1
        end
         
        # After the first `while` is completed, we need to change the value of  the modifying constant `f`.
        lnf= (lnf/2) # Update `f`.
    end
    
    return enDensityDict
end















"""
    determineλ(lngE,T)
Given a dictionary `lngE` whose keys(values) are energies(log of energy densities), and a temperature `T`; returns the largest value
λ of the differences between energy density `ln(g(E))` and `βE` (i.e. `max{log(g(Eᵢ))-βE}`).
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
    sExp(E,T,lnge,λ)
Given a temperature `T`, an energy `E`, the logarithm of the enegy density `lnge` and a number `λ`; returns e^{ln(g(e))-βE-λ}.
"""
function sExp(E,T,lnge,λ)
    β = 1/T
    res = exp(lnge-β*E-λ)
    return res
end










"""
    main_WL(N,protein,ti,tf,nTs,numlim2,pfmodel::PF_model,name)
 Given a lattice size `N`, a protein `protein`, an initial (final) temperature `ti` (`tf`), the number of visited temperatures `nTs`, 
 a limit number of iterations `numlim2` and an interaction model `pfmodel`, and a name where to store the data `name`; stores the 
 values of the thermodynamic variables to the desired file.   
"""
function main_WL(N,protein,ti,tf,nTs,numlim2,pfmodel::PF_model,name::String)

    edo = protein.edo
    HPlist = protein.HPlist
    geometry = protein.geometry
    temps = range(ti,stop=tf,length=nTs)

    # Create the directory which will contain the dat collected trough the simulation.
    pathstring="./outputWL/"
    pathname=pathstring*name
    mkdir(pathname)
    writedlm(pathname*"/temperatures.csv",temps,',')
    writedlm(pathname*"/initialconf.csv",edo,',') # Perhaps a tad ineccesary.
    writedlm(pathname*"/HPlist.csv",Int.(HPlist),',')
    writedlm(pathname*"/latticesize.csv",N,',') # Perhaps a tad ineccesary.
    writedlm(pathname*"/geometry.csv",Int(protein.geometry),',') # Need to import the translation function.
    writedlm(pathname*"/pfmodel.csv",Int(pfmodel.pf_name),',') # Need to import the dictionary that interprets.
 
    us = ones(Float64,nTs) # This array will contain the internal energy per aminoacid for each of the temperatures in `temps`.
    cs=ones(Float64,nTs) 
    Fs = ones(Float64,nTs)
    Ss = ones(Float64,nTs)

    lngE = wang_landau(N,edo,HPlist,geometry,numlim2,pfmodel) # Compute the natural logarithm of the energy densities.
    
    # I need to normalize the density of states.
    #### !!!!! Need help here !!!!!!!!!!!!
    # mink = minimum(keys(lngEaux)) # Maximum energy recorded during the simulation.
    # maxk = maximum(keys(lngEaux)) # Minimum energy recorded during the simulation.
    # mrest = min(lngEaux[mink],lngEaux[maxk]) # Choose the minimum density of states.
    mrest = minimum(Float64[val for val in values(lngE)]) # Choose the minimum density of states.
    
    # Declare a normalized dictionary.
    lngE = Dict{Float64,Float64}(k => (lngE[k]-mrest) for k in keys(lngE))

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
            𝚭 = 𝚭 + z
        end
        
        uT = (𝚬/𝚭) # Internal energy for the given temperature.
        cT = ((𝚬sq/𝚭)-(uT^2))/(T^2) 
        fT = -T*(λ+log(𝚭))
        entropyS = (uT-fT)/T
        
        us[i] = uT/N # Energy per spin.
        cs[i] = cT/N # Specific heat per spin.
        Fs[i] = fT/N
        Ss[i] = entropyS/N
    end


    writedlm(pathname*"/us.csv",us,',')
    writedlm(pathname*"/cs.csv",cs,',')
    writedlm(pathname*"/Fs.csv",Fs,',')
    writedlm(pathname*"/Ss.csv",Ss,',')
    writedlm(pathname*"/lngE.csv",lngE,',')

    println("Wang-Landau simulation succesfully completed.")
end


