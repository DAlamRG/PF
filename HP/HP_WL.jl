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
function energyDif(N::Int,edoi::Array{Int16,2},edof::Array{Int16,2},HPlist::Vector{Amin},geometry::geometries,pfmodel::PF_model)
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
    wang_landau(N,protein,numlim2,pfmodel::PF_model,name)
Given a lattice size `N`, a protein structure `protein`, a limit number of iterations `numlim2`, an interaction model `pfmodel` and
a name where for the fle where the dta will be stored `name`; returns a dictionary whose keys are the visited energies, and 
it's values are the logarithms of the energy densities.
"""
function wang_landau(N::Int,protein::Protein,numlim2::Int,pfmodel::PF_model,name::String)
    
    edo = protein.edo
    HPlist = protein.HPlist
    geometry = protein.geometry

    mcsweep = length(edo) # Monte Carlo step size.
    
    lnf = 1.0 # Initial value of updating value `log(f)`.
    enDensityDict = Dict{Float64,Float64}(0 => 0) # Dictionary for the energy densities. Keys are all the possible energy values
    # for the system (I start with just one since the range of possible energies depends entirely on the aminoacid list/geometry/interaction model).
    # Values are the energy density (actually log(g(Eᵢ)) ), where g(Eᵢ) is the energy density.


    # Before performing the actual simulation, I'd like to find most of the possible energies. This is crucial.

    es = Float64[]
    for k in 1:(5000*mcsweep) # Perform 5000 Monte Carlo steps.
        e1 = energy(N,edo,HPlist,geometry,pfmodel) # Initial energy.
        npull_1 = countpull(N,edo,HPlist,geometry)[1]
        # Generate new state by performing a pull-move.
        edoaux,npull_2 = pullMove(N,edo,HPlist,geometry)[1:2]
        e2 = energy(N,edoaux,HPlist,geometry,pfmodel) # Energy of the new state

        # There´s a chance that the computed energies are not yet in the dictionary. In that case, we add them.
        if e1 ∉ keys(enDensityDict)
            enDensityDict[e1] = 0
        end
        if e2 ∉ keys(enDensityDict)
            enDensityDict[e2] = 0
        end

        # Next, we need to obtain the density of states for the energies.
        lg1 = enDensityDict[e1]
        lg2 = enDensityDict[e2]

        r = rand() # Generate a random number.
        pμν = exp(lg1-lg2+log(npull_1)-log(npull_2))
        if r ≤ pμν # Accept the flip.
            edo = edoaux # Accept the new configuration.  
        end            

    end

    for e in minimum(keys(enDensityDict)):0
        if e ∉ keys(enDensityDict)
            enDensityDict[e2] = 0
        end

    end

    @show(enDensityDict)




    # Now I do perform the actual simulation. 

    cont1 = 1
    for l in 1:27 # This is the number of iterations it takes to make lnf sufficiently small.
        println("cont1= $cont1 /27")
        cont1 = cont1+1 # Update the counter.
        localenergies = Dict{Float64,Int64}(k => 0 for k in keys(enDensityDict)) # Stores the number of times each energy is visted during the current iteration. This is my "histogram".
        
        hCond = false # True if the histogram is flat enough.
        cont2 = 1
        while (cont2 ≤ numlim2) && (hCond == false)

            for k in 1:(100*mcsweep) # Perform 50 Monte Carlo steps.
                e1 = energy(N,edo,HPlist,geometry,pfmodel) # Initial energy.
                npull_1 = countpull(N,edo,HPlist,geometry)[1]
                # Generate new state by performing a pull-move.
                edoaux,npull_2 = pullMove(N,edo,HPlist,geometry)[1:2]
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
                pμν = exp(lg1-lg2+log(npull_1)-log(npull_2))
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

    # Next, we need to normalize the dictionary.
    mrest = minimum(Float64[val for val in values(enDensityDict)]) # Choose the minimum density of states.
    
    # Declare a normalized dictionary.
    lngE = Dict{Float64,Float64}(k => (enDensityDict[k]-mrest) for k in keys(enDensityDict))

    # Create the directory which will contain the data collected trough the simulation.
    pathstring = "./outputWL/"
    pathname = pathstring*name
    mkdir(pathname)
    writedlm(pathname*"/initialconf.csv",edo,',') # Perhaps a tad ineccesary.
    writedlm(pathname*"/HPlist.csv",Int.(HPlist),',')
    writedlm(pathname*"/latticesize.csv",N,',') 
    writedlm(pathname*"/geometry.csv",Int(protein.geometry),',') # Need to import the translation function.
    writedlm(pathname*"/pfmodel.csv",Int(pfmodel.pf_name),',') # Need to import the dictionary that interprets.
    writedlm(pathname*"/lngE.csv",lngE,',')
    
    println("WL simulation completed.")
end

