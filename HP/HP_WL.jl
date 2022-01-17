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
    histCondtion(localenergies,greek)

Given a dictionary `localenergies` containing the number of times each energy was visited and a percentage `percent`; counts the incidence of each energy. If the
lowest counted energy is `percent` % of the average, return a boolean value `true`.
"""
function histCondtion(localenergies,percent,greek)
    
    if greek ==  true
        ocurrences = values(localenergies)
        hmin = minimum(ocurrences)
        m = mean(ocurrences)
        err = (m-hmin)/m
        val = hmin ≥ m*(percent*1e-2)
        return (val,err)
    else
        ocurrences = localenergies
        hmin = minimum(ocurrences)
        m = mean(ocurrences)
        err = (m-hmin)/m
        val = hmin ≥ m*(percent*1e-2)
        return (val,err)
    end
    
end













"""
    metropolis(N,nums,T,protein,pfmodel) 

Given a 2D/3D array size `N`, a number of Monte-Carlo sweeps `nums`, the temperature T, a protein structure `protein` encoding the protein's 
configuration and an amino acid interaction model `pfmodel`; returns the information necessay to reproduce the visited configuration during the 
simulation (using the Metropolis-Hastings algorithm) as well as the visited energies.
"""
function metropolis(N::Int,nums::Int,T::Float64,protein::Protein,pfmodel::PF_model)  
    
    β = 1/T
    edo = protein.edo
    HPlist = protein.HPlist
    geometry = protein.geometry
    ns = nums*length(HPlist) # Total number of iterations

    prev_state = edo # This variable stores the previous configuration.
    prev_numpull = countpull(N,edo,HPlist,geometry)[1] # This variable stores the number of possible pull moves for the previous iteration.
    
    enstates = zeros(Float64,ns+1) # Stores the energy of each state visited during the simulation.
    enstates[1] = energy(N,edo,HPlist,geometry,pfmodel) # Compute the energy of the initial state and store it.
    
    for l in 2:ns+1 
        # Generate a new state.
        newedo,totalpull,indpulled,newcoord,dirpulled = pullMove(N,prev_state,HPlist,geometry)
        newenergy = energy(N,newedo,HPlist,geometry,pfmodel) # Compute the energy after the pull move.
        ΔH = newenergy-enstates[l-1] # Compute the energy difference.
        
        if ΔH ≤ 0 # We accept the configuration change.
            prev_state = newedo
            enstates[l] = newenergy 
            prev_numpull = totalpull

        else
            r = rand() # Random number, sampled from uniform distribution over [0,1].
            exponential = exp(-β*ΔH)
            q = (prev_numpull)/(totalpull) # This is the quotient of number of pull moves.
            pστ = q*exponential # This determines the probability of the move in the case ΔH > 0.

            if r < pστ # Accept the move.
                prev_state = newedo
                enstates[l] = newenergy 
                prev_numpull = totalpull

            else # Move is not accepted, the protein´s configuration stays the same.
                enstates[l] = copy(enstates[l-1])
                
            end
        end
    end

    return (enstates,prev_state)

end












"""
    prev_met(N,nums,ti,tf,nTs,protein,pfmodel)

Given a 2D/3D array size `N`, a number of Monte-Carlo sweeps `nums`, the initial/final temperature `ti(tf)`, a number of temperatures `nTs`,
a protein structure `protein`, and an amino acid interaction model `pfmodel`; returns an array with all of the visited energies.
"""
function prev_met(N::Int,nums::Int,ti::Float64,tf::Float64,nTs::Int,protein::Protein,pfmodel::PF_model)
    temperatures = range(ti,stop=tf,length=nTs) # Declare a range of temperatures to be visited.
    ns = nums*length(protein.HPlist)
    energies = zeros(Float64,Int((ns*length(temperatures))+1))

    # Perform the first simulation.
    enstates,laststate = metropolis(N,nums,temperatures[end],protein,pfmodel)
    energies[1:ns+1] = enstates
    
    # For each of the remaining temperatures, employ the Metropolis-Hastings algorithm to store the information about the
    # visited configurations at the current temperature.
    cont = ns+2
    for k in 2:length(temperatures)
        temp = reverse(temperatures)[k] # temperature at which the simulation is performed
        proteinaux = Protein(laststate,protein.HPlist,protein.geometry) # Protein structure for the simulation.
        enstates,laststate = metropolis(N,nums,temp,proteinaux,pfmodel)  # Perfom the simulation.
        energies[cont:cont+(ns-1)] = enstates[2:end] # Store the visited energies.
        cont = cont+ns 
    end

    return energies
        
end









"""
    energy_bins(e_min,e_max,pfmodel)
Given a value for the lowest/greatest visited energy `e_min/e_max` and a folding model `pfmodel`; returns binned energies and an array for the density of states.
"""
function energy_bins(e_min,e_max,pfmodel::PF_model)
    pf_M = pfmodel.int_matrix
    aux_vec = setdiff!(abs.(vcat(pf_M...)), [0])
    step_size = maximum(aux_vec)
    m = abs(e_min)+(0.0000001)+abs(e_max)
    n_bins = Int(ceil(m/step_size))
    energies = reverse(collect(0.0000001+abs(e_max):-step_size:-n_bins*step_size)) 
    return (energies,zeros(Float64,length(energies)-1))
end







"""
     update_energy_bins!(e_new,num,energies,enDensityDict,pfmodel)
Given a new minimum/maximum energy `e_new`, a value `num` (a value of `1/2` indicates that the new energy is lower/greater), the previous binned energies `energies`, the array containing the density of states `enDensityDict`,
and a folding model; returns updated binned energies.
"""
function update_energy_bins!(e_new,num,energies,enDensityDict,pfmodel::PF_model)
    pf_M = pfmodel.int_matrix
    aux_vec = setdiff!(abs.(vcat(pf_M...)), [0])
    step_size = maximum(aux_vec)

    if num == 1
        m = abs(e_new)+minimum(energies)
        n_bins = Int(ceil(m/step_size))
        energies1 = collect(range(-n_bins*step_size + minimum(energies),stop=minimum(energies),step = step_size))[1:end-1]
        energies2 = vcat(energies1,energies)
        enDensityDict2 = vcat(fill(minimum(enDensityDict),length(energies1)),enDensityDict)
        return (energies2,enDensityDict2)
    else
        m = abs(e_new)-maximum(energies)
        n_bins = Int(ceil(m/step_size))
        energies1 = collect(range(maximum(energies),stop=maximum(energies)+n_bins*step_size,step = step_size))[2:end]
        energies2 = vcat(energies,energies1)
        enDensityDict2 = vcat(enDensityDict,fill(minimum(enDensityDict),length(energies1)))
        return (energies2,enDensityDict2)
    end
end









#=
"""
    wang_landau_beta(N,protein,numlim2,pfmodel::PF_model,name)
Given a lattice size `N`, a protein structure `protein`, a limit number of iterations `numlim2`, an interaction model `pfmodel` and
a name where for the fle where the dta will be stored `name`; returns a dictionary whose keys are the visited energies, and 
it's values are the logarithms of the energy densities.
"""
function wang_landau_beta(N::Int,protein::Protein,numlim2::Int,pfmodel::PF_model,name::String)
    
    edo = protein.edo
    HPlist = protein.HPlist
    geometry = protein.geometry

    mcsweep = length(edo) # Monte Carlo step size.
    
    lnf = 1.0 # Initial value of updating value `log(f)`.

    # Before performing the actual simulation, I'd like to find most of the possible energies. This is crucial.
    # To do this, I perform a Metropolis simulation over a range of adequate temperatures.

    es = prev_met(N,100,0.01,2.5,20,protein,pfmodel)
    es = unique(es)
    @show(es)

    enDensityDict = Dict{Float64,Float64}(e => 0 for e in es) # Dictionary for the energy densities. Keys are all the possible energy values
    # for the system (I start with just one since the range of possible energies depends entirely on the aminoacid list/geometry/interaction model).
    # Values are the energy density (actually log(g(Eᵢ)) ), where g(Eᵢ) is the energy density.


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
                e1 = round(energy(N,edo,HPlist,geometry,pfmodel), digits = 3) # Initial energy.
                npull_1 = countpull(N,edo,HPlist,geometry)[1]
                # Generate new state by performing a pull-move.
                edoaux,npull_2 = pullMove(N,edo,HPlist,geometry)[1:2]
                e2 = round(energy(N,edoaux,HPlist,geometry,pfmodel), digits = 3) # Energy of the new state

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
    
                lr = log(rand()) # Generate a random number.
                pμν = lg1-lg2+log(npull_1)-log(npull_2)
                if lr ≤ pμν # Accept the flip.
                    edo = edoaux # Accept the new configuration.  
                    enDensityDict[e2] = lg2+lnf # Update the density of states for the new energy (the log actually).
                    localenergies[e2] = localenergies[e2] + 1 # Update the "histogram".
    
                else
                    enDensityDict[e1] = lg1+lnf # Update the density of states for the new energy (the log actually).
                    localenergies[e1] = localenergies[e1] + 1  # Update the "histogram".
                end            
    
            end

            hCond = histCondtion(localenergies,true)
            if hCond == true 
                println("hCond = ",hCond)
            end

            cont2 = cont2 + 1
        end
         
        # After the first `while` is completed, we need to change the value of  the modifying constant `f`.
        lnf = (lnf/2) # Update `f`.
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
=#














"""
    wang_landau(N,protein,numlim2,pfmodel,name, energy_bins)
Given a lattice size `N`, a protein structure `protein`, a limit number of iterations `numlim2`, an interaction model `pfmodel`,
a name where for the fle where the data will be stored `name`; returns a dictionary whose 
keys are the visited energies, and it's values are the logarithms of the energy densities.
"""
function wang_landau(N::Int,protein::Protein,numlim2::Int,pfmodel::PF_model,name::String)
    
    edo = protein.edo
    HPlist = protein.HPlist
    geometry = protein.geometry

    mcsweep = length(edo) # Monte Carlo step size.
    
    lnf = 1.0 # Initial value of updating value `log(f)`.

    # Before performing the actual simulation, I'd like to find most of the possible energies. This is crucial.
    # To do this, I perform a Metropolis simulation over a range of adequate temperatures.
    es = prev_met(N,170,0.001,2.5,24,protein,pfmodel)
    unique!(es)
    sort!(es)
    @show(es)
    energies, enDensityDict = energy_bins(minimum(es),maximum(es),pfmodel)
    
     
    # Now I perform the actual simulation. 

    cont1 = 1
    hCond_vec = fill(false,27)
    err_vec = zeros(Float64,27)
    for l in 1:27 # This is the number of iterations it takes to make lnf sufficiently small.
        println("cont1= $cont1 /27")
        cont1 = cont1+1 # Update the counter.
        localenergies = zeros(Int64,length(enDensityDict)) # Stores the number of times each energy is visted during the current iteration.
        
        hCond = false # True if the histogram is flat enough.
        err = 10.0
        cont2 = 1
        while (cont2 ≤ numlim2) && (hCond == false)

            for k in 1:(1000*mcsweep) # Perform 1000 Monte Carlo steps.
                e1 = energy(N,edo,HPlist,geometry,pfmodel) # Initial energy.
                npull_1 = countpull(N,edo,HPlist,geometry)[1]
                # Generate new state by performing a pull-move.
                edoaux,npull_2 = pullMove(N,edo,HPlist,geometry)[1:2]
                e2 = energy(N,edoaux,HPlist,geometry,pfmodel) # Energy of the new state

                if e2 < energies[1]
                    energies, enDensityDict = update_energy_bins!(e2,1,energies,enDensityDict,pfmodel)
                    @show(energies)
                    @show(enDensityDict)
                    localenergies = zeros(Int64,length(enDensityDict))
                    cont2 = 1
                elseif e2 > energies[end]
                    energies, enDensityDict = update_energy_bins!(e2,2,energies,enDensityDict,pfmodel)
                    @show(energies)
                    @show(enDensityDict)
                    localenergies = zeros(Int64,length(enDensityDict))
                    cont2 = 1
                end

                # Next, we need to obtain the density of states for the energies.
                ind_e1 = searchsortedlast(energies,e1)
                ind_e2 = searchsortedlast(energies,e2)
                lg1 = enDensityDict[ind_e1]
                lg2 = enDensityDict[ind_e2]
    
    
                lr = log(rand()) # Generate a random number.
                pμν = lg1-lg2+log(npull_1)-log(npull_2)
                if lr ≤ pμν # Accept the flip.
                    edo = edoaux # Accept the new configuration.  
                    enDensityDict[ind_e2] = lg2+lnf # Update the density of states for the new energy (the log actually).
                    localenergies[ind_e2] = localenergies[ind_e2]+1 # Update the "histogram".
    
                else
                    enDensityDict[ind_e1] = lg1+lnf # Update the density of states for the new energy (the log actually).
                    localenergies[ind_e1] = localenergies[ind_e1]+1  # Update the "histogram".
                end            
    
            end

            hCond, err = histCondtion(localenergies,80,false)
            if hCond == true 
                println("hCond = ",hCond)
            end

            cont2 = cont2 + 1
        end
         
        # After the first `while` is completed, we need to change the value of  the modifying constant `f`.
        lnf = (lnf/2) # Update `f`.
        hCond_vec[l] = hCond
        err_vec[l] = err
    end

    @show(enDensityDict)

    # Next, we need to normalize the dictionary.
    mrest = minimum(enDensityDict) # Choose the minimum density of states.
    
    # Declare a normalized dictionary.
    lngE = Dict{Float64,Float64}(energies[k+1] => (enDensityDict[k]-mrest) for k in 1:length(energies)-1)
    # Dict{Float64,Float64}((energies[k]+energies[k+1])/2 => (enDensityDict[k]-mrest) for k in 1:length(energies)-1)
    @show(lngE)
    # Create the directory which will contain the data collected trough the simulation.
    pathstring = "./outputWL/"
    pathname = pathstring*name
    mkdir(pathname)
    writedlm(pathname*"/initialconf.csv",edo,',') # Perhaps a tad ineccesary.
    writedlm(pathname*"/HPlist.csv",Int.(HPlist),',')
    writedlm(pathname*"/latticesize.csv",N,',')
    writedlm(pathname*"/hCond_vec.csv",Int.(hCond_vec),',')
    writedlm(pathname*"/err_vec.csv",err_vec,',')  
    writedlm(pathname*"/numlim2.csv",numlim2,',') 
    writedlm(pathname*"/geometry.csv",Int(protein.geometry),',') # Need to import the translation function.
    writedlm(pathname*"/pfmodel.csv",Int(pfmodel.pf_name),',') # Need to import the dictionary that interprets.
    writedlm(pathname*"/lngE.csv",lngE,',')
    
    println("WL simulation completed.")
end









##### This alternative functions are to be used if I want to use the original interaction matrices in the simulations.
"""
    energy_bins(energies,pfmodel)
Given a vector of the visited energy levels `energies` and a folding model `pfmodel`; returns binned energies and an array for the density of states.
"""
function energy_bins(energies,pfmodel::PF_model)
    pf_n = pfmodel.pf_name

    if pf_n == Full1 || pf_n == Full2
        e_min = min(energies) 
        e_max = max(energies)
        pf_M = pfmodel.int_matrix
        aux_vec = setdiff!(abs.(vcat(pf_M...)), [0])
        step_size = maximum(aux_vec)
        m = abs(e_min)+(0.0000001)+abs(e_max)
        n_bins = Int(ceil(m/step_size))
        energies = reverse(collect(0.0000001+abs(e_max):-step_size:-n_bins*step_size)) 
        return (energies,zeros(Float64,length(energies)-1))

    else
        unique!(energies)
        sort!(energies)
        return (energies,zeros(Float64,length(energies)))

    end
end



"""
     update_energy_bins!(e_new,num,energies,enDensityDict,pfmodel)
Given a new minimum/maximum energy `e_new`, a value `num` (a value of `1/2` indicates that the new energy is lower/greater), the previous binned energies `energies`, the array containing the density of states `enDensityDict`,
and a folding model; returns updated binned energies.
"""
function update_energy_bins!(e_new,num,energies,enDensityDict,pfmodel::PF_model)
    pf_n = pfmodel.pf_name

    if pf_n == Full1 || pf_n == Full2 
        
        pf_M = pfmodel.int_matrix
        aux_vec = setdiff!(abs.(vcat(pf_M...)), [0])
        step_size = maximum(aux_vec)
        if num == 1
            m = abs(e_new)+minimum(energies)
            n_bins = Int(ceil(m/step_size))
            energies1 = collect(range(-n_bins*step_size + minimum(energies),stop=minimum(energies),step = step_size))[1:end-1]
            energies2 = vcat(energies1,energies)
            enDensityDict2 = vcat(fill(minimum(enDensityDict),length(energies1)),enDensityDict)
            return (energies2,enDensityDict2)
        else
            m = abs(e_new)-maximum(energies)
            n_bins = Int(ceil(m/step_size))
            energies1 = collect(range(maximum(energies),stop=maximum(energies)+n_bins*step_size,step = step_size))[2:end]
            energies2 = vcat(energies,energies1)
            enDensityDict2 = vcat(enDensityDict,fill(minimum(enDensityDict),length(energies1)))
            return (energies2,enDensityDict2)
        end

    else

        if num == 1
            pushfirst!(energies,e_new)
            enDensityDict2 = vcat(fill(minimum(enDensityDict),1),enDensityDict)
            return (energies,enDensityDict2)
        else
            push!(energies,e_new)
            enDensityDict2 = vcat(enDensityDict,fill(minimum(enDensityDict),1))
            return (energies,enDensityDict2)
        end
        
    end
end