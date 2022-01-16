using Base: Int16, Float64, String
# This script performs a Metropolis simulation, and stores the data.
using Distributed 
Distributed.addprocs(4)



@sync @everywhere begin
using DelimitedFiles

include("./Energy.jl")



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
    
    
    # Now, apply Metropolis-Hastings.
    
    if geometry == fcc || geometry == triangular2D
        pulledindices = zeros(Int16,ns) # Stores the moved indices. If no index is moved, the value is zero.
        newcoords = zeros(Int16,(ns,length(edo[1,:]))) # Stores the new coordinate for the pulled index. If the configuration stays the same, the new coordinate is just zeros.
        dirs = fill(nonetaken,ns) # Store the direction in which the chain was pulled. If no pull move was performed, the direction is `nonetaken`.

        for l in 2:ns+1 

            # Generate a new state.
            newedo,totalpull,indpulled,newcoord,dirpulled = pullMove(N,prev_state,HPlist,geometry)
            newenergy = energy(N,newedo,HPlist,geometry,pfmodel) # Compute the energy after the pull move.
            
            ΔH = newenergy-enstates[l-1] # Compute the energy difference.
            
    
            if ΔH ≤ 0 # We accept the configuration change.
                prev_state = newedo
                enstates[l] = newenergy 
                prev_numpull = totalpull
                pulledindices[l-1] = indpulled
                newcoords[l-1,:] = newcoord
                dirs[l-1] = dirpulled
    
            else
                r = rand() # Random number, sampled from uniform distribution over [0,1].
                exponential = exp(-β*ΔH)
                q = (prev_numpull)/(totalpull) # This is the quotient of number of pull moves.
                pστ = q*exponential # This determines the probability of the move in the case ΔH > 0.
    
                if r < pστ # Accept the move.
                    prev_state = newedo
                    enstates[l] = newenergy 
                    prev_numpull = totalpull
                    pulledindices[l-1] = indpulled
                    newcoords[l-1,:] = newcoord
                    dirs[l-1] = dirpulled
    
                else # Move is not accepted, the protein´s configuration stays the same.
                    enstates[l] = copy(enstates[l-1])
                    
                end
            end
        end

        return (pulledindices,dirs,newcoords,enstates)
        

    elseif geometry == square2D
        pulledindices = zeros(Int16,ns) # Stores the moved indices. If no index is moved, the value is zero.
        newcoords1 = zeros(Int16,(ns,length(edo[1,:]))) # Stores the new coordinate for the pulled index. 
        newcoords2 = zeros(Int16,(ns,length(edo[1,:]))) # Stores the new coordinate for the momomer directly before/after the one being moved. 
        dirs = fill(nonetaken,ns) # Store the direction in which the chain was pulled. If no pull move was performed, the direction is `nonetaken`.

        #Generate new states by performing pull moves.
        for l in 2:ns+1 
            # Generate a new state.
            newedo,totalpull,indpulled,newcoord1,newcoord2,dirpulled = pullMove(N,prev_state,HPlist,geometry)
            newenergy = energy(N,newedo,HPlist,geometry,pfmodel) # Compute the energy after the pull move.
            
            ΔH = newenergy-enstates[l-1] # Compute the energy difference.
            
    
            if ΔH ≤ 0 # We accept the configuration change.
                prev_state = newedo
                enstates[l] = newenergy 
                prev_numpull = totalpull
                pulledindices[l-1] = indpulled
                newcoords1[l-1,:] = newcoord1
                newcoords2[l-1,:] = newcoord2
                dirs[l-1] = dirpulled
    
    
            else
                r = rand() # Random number, sampled from uniform distribution over [0,1].
                exponential = exp(-β*ΔH)
                q = (prev_numpull)/(totalpull) # This is the quotient of number of pull moves.
                pστ = q*exponential # This determines the probability of the move in the case ΔH > 0.
    
                if r < pστ # Accept the move.
                    prev_state = newedo
                    enstates[l] = newenergy 
                    prev_numpull = totalpull
                    pulledindices[l-1] = indpulled
                    newcoords1[l-1,:] = newcoord1
                    newcoords2[l-1,:] = newcoord2
                    dirs[l-1] = dirpulled

                else # Move is not accepted, the protein´s configuration stays the same.
                    enstates[l] = copy(enstates[l-1])

                end
            end
        end

        newcoords = (newcoords1,newcoords2)
        return (pulledindices,dirs,newcoords,enstates)

    end

end

end



















 # Next, I write a function which performs multiple simulations over an array of temperatures, and saves the generated info.
"""
     main_met_1(N,nums,ti,tf,nTs,nruns,protein,pfmodel,name)

Given a 2D/3D array size `N`, a number of Monte-Carlo sweeps per temperature `nums`, an initial(final) temperature `ti(tf)`, the number of 
temperatures to be visited `nTs`, a number of independent temperature sweeps `nruns`, a structure `protein` encoding the protein´s 
configuration, an aminoacid interaction model `pfmodel` and a name for the simulation output `name`; writes the information needed to 
recreate the visited states after performing a simulation using the Metropolis-Hastings algorithm.
"""
function main_met_1(N::Int,nums::Int,ti::Float64,tf::Float64,nTs::Int,nruns::Int,protein::Protein,pfmodel::PF_model,name::String)
    temperatures = range(ti,stop=tf,length=nTs) # Declare a range of temperatures to be visited.
    ns = nums*length(protein.HPlist)

    # Create the directory which will contain the dat collected trough the simulation.
    pathstring = "./output/"
    pathname = pathstring*name
    mkdir(pathname)
    writedlm(pathname*"/temperatures.csv",temperatures,',')
    writedlm(pathname*"/initialconf.csv",protein.edo,',')
    writedlm(pathname*"/HPlist.csv",Int.(protein.HPlist),',')
    writedlm(pathname*"/latticesize.csv",N,',')
    writedlm(pathname*"/mc_sweeps.csv",nums,',')
    writedlm(pathname*"/geometry.csv",Int(protein.geometry),',')
    writedlm(pathname*"/pfmodel.csv",Int(pfmodel.pf_name),',')
    writedlm(pathname*"/nruns.csv",nruns,',')

    geometry = protein.geometry

    if geometry == fcc || geometry == triangular2D
        # Now, perform a Metropolis-Hastings simulation for each temperature in `temperatures`. Sweep the tempeartures `nruns` times.
        for l in 1:nruns
            energies = zeros(Float64,Int((ns*length(temperatures))+1)) #Float64[] # Energies for the current run.
        
            # Perform the first simulation.
            datatemp1 = metropolis(N,nums,temperatures[end],protein,pfmodel)
            pulledindicesT = datatemp1[1] # The next three variables will be used to perform subsequent simulations.
            dirsT = datatemp1[2]
            newcoordsT = datatemp1[3]
            laststate = reconstructStates(N,protein.edo,protein.HPlist,pulledindicesT,dirsT,newcoordsT,protein.geometry)[:,:,end]
            energies[1:ns+1] = datatemp1[4]
            # append!(energies,datatemp1[4]) # Store the first batch of energies.
        
        
            # Store the output generated in the first simulation.
            pathnameaux = pathname*"/"*string(l)
            writedlm(pathnameaux*"_1_1.csv",pulledindicesT,',')
            writedlm(pathnameaux*"_1_2.csv",Int.(dirsT)) # Save only the numbers, not the whole enum type.
            writedlm(pathnameaux*"_1_3.csv",newcoordsT,',')
        
            # For each of the remaining temperatures, employ the Metropolis-Hastings algorithm to store the information about the
            # visited configurations at the current temperature.
            cont = ns+2
            for k in 2:length(temperatures)
                temp = reverse(temperatures)[k] # temperature at which the simulation is performed
                proteinaux = Protein(laststate,protein.HPlist,protein.geometry) # Protein structure for the simulation.
                pulledindicesT,dirsT,newcoordsT,enstates = metropolis(N,nums,temp,proteinaux,pfmodel)  # Perfom the simulation.
                # append!(energies,enstates) # Store the visited energies.
                energies[cont:cont+(ns-1)] = enstates[2:end] # Store the visited energies.
                cont = cont+ns 
                
                
                # Store the output generated in the first simulation.
                st = "_"*string(k)
                writedlm(pathnameaux*st*"_1.csv",pulledindicesT,',')
                writedlm(pathnameaux*st*"_2.csv",Int.(dirsT),',')
                writedlm(pathnameaux*st*"_3.csv",newcoordsT,',')
                
                laststate = reconstructStates(N,laststate,protein.HPlist,pulledindicesT,dirsT,newcoordsT,protein.geometry)[:,:,end] 
                println("Temperature Progress: $k /$nTs")
            end
            writedlm(pathnameaux*"_energies.csv",energies,',') # Save all of the visted energies.
            println("Number of runs progress : $l /$nruns", )
        end
        
        
    elseif geometry == square2D
        # Now, perform a Metropolis-Hastings simulation for each temperature in `temperatures`. Sweep the tempeartures `nruns` times.
        for l in 1:nruns
            energies = zeros(Float64,Int((ns*length(temperatures))+1)) # Energies for the current run.

            # Perform the first simulation.
            datatemp1 = metropolis(N,nums,temperatures[end],protein,pfmodel)
            pulledindicesT = datatemp1[1] # The next three variables will be used to perform subsequent simulations.
            dirsT = datatemp1[2]
            newcoordsT = datatemp1[3]
            laststate = reconstructStates(N,protein.edo,protein.HPlist,pulledindicesT,dirsT,newcoordsT,protein.geometry)[:,:,end]
            energies[1:ns+1] = datatemp1[4] # Store the first batch of energies.
    
            # Store the output generated in the first simulation.
            pathnameaux = pathname*"/"*string(l)
            writedlm(pathnameaux*"_1_1.csv",pulledindicesT,',')
            writedlm(pathnameaux*"_1_2.csv",Int.(dirsT),',') # Save only the Int value, not the whole enum info.
            writedlm(pathnameaux*"_1_3_1.csv",newcoordsT[1],',')
            writedlm(pathnameaux*"_1_3_2.csv",newcoordsT[2],',')

    
            # For each of the remaining temperatures, employ the Metropolis-Hastings algorithm to store the information about the
            # visited configurations at the current temperature.
            cont = ns+2
            for k in 2:length(temperatures)
                temp = reverse(temperatures)[k] # temperature at which the simulation is performed
                proteinaux = Protein(laststate,protein.HPlist,protein.geometry) # Auxiliary protein structure for the simulation.
                pulledindicesT,dirsT,newcoordsT,enstates = metropolis(N,nums,temp,proteinaux,pfmodel)  # Perform the simulation.
                energies[cont:cont+(ns-1)] = enstates[2:end] # Store the visited energies.
                cont = cont+ns 

                # Store the output generated in the first simulation.
                st = "_"*string(k)
                writedlm(pathnameaux*st*"_1.csv",pulledindicesT,',')
                writedlm(pathnameaux*st*"_2.csv",Int.(dirsT),',')
                writedlm(pathnameaux*st*"_3_1.csv",newcoordsT[1],',')
                writedlm(pathnameaux*st*"_3_2.csv",newcoordsT[2],',')
                
                laststate = reconstructStates(N,laststate,protein.HPlist,pulledindicesT,dirsT,newcoordsT,protein.geometry)[:,:,end] 
                println("Temperature Progress: $k /$nTs")
            end

            writedlm(pathnameaux*"_energies.csv",energies,',') # Save all of the visted energies.
            println("Number of runs progress : $l /$nruns", )

        end

    end
    
    println("Simulation succesfully completed")
end










"""
     main_met(N,nums,ti,tf,nTs,nruns,protein,pfmodel,name)

Given a 2D/3D array size `N`, a number of Monte-Carlo sweeps per temperature `nums`, an initial(final) temperature `ti(tf)`, the number of 
temperatures to be visited `nTs`, a number of independent temperature sweeps `nruns`, a structure `protein` encoding the protein´s 
configuration, an aminoacid interaction model `pfmodel` and a name for the simulation output `name`; writes the information needed to 
recreate the visited states after performing a simulation using the Metropolis-Hastings algorithm.
"""
function main_met(N::Int,nums::Int,ti::Float64,tf::Float64,nTs::Int,nruns::Int,protein::Protein,pfmodel::PF_model,name::String)
    temperatures = range(ti,stop=tf,length=nTs) # Declare a range of temperatures to be visited.
    ns = nums*length(protein.HPlist)

    # Create the directory which will contain the dat collected trough the simulation.
    pathstring = "./output/"
    pathname = pathstring*name
    mkdir(pathname)
    mkdir(pathname*"/energiesFolder")
    writedlm(pathname*"/temperatures.csv",temperatures,',')
    writedlm(pathname*"/initialconf.csv",protein.edo,',')
    writedlm(pathname*"/HPlist.csv",Int.(protein.HPlist),',')
    writedlm(pathname*"/latticesize.csv",N,',')
    writedlm(pathname*"/mc_sweeps.csv",nums,',')
    writedlm(pathname*"/geometry.csv",Int(protein.geometry),',')
    writedlm(pathname*"/pfmodel.csv",Int(pfmodel.pf_name),',')
    writedlm(pathname*"/nruns.csv",nruns,',')

    geometry = protein.geometry

    if geometry == fcc || geometry == triangular2D
        # Now, perform a Metropolis-Hastings simulation for each temperature in `temperatures`. Sweep the tempeartures `nruns` times.
        @sync @distributed for l in 1:nruns
            energies = zeros(Float64,Int((ns*length(temperatures))+1)) #Float64[] # Energies for the current run.
        
            # Perform the first simulation.
            datatemp1 = metropolis(N,nums,temperatures[end],protein,pfmodel)
            pulledindicesT = datatemp1[1] # The next three variables will be used to perform subsequent simulations.
            dirsT = datatemp1[2]
            newcoordsT = datatemp1[3]
            laststate = reconstructStates(N,protein.edo,protein.HPlist,pulledindicesT,dirsT,newcoordsT,protein.geometry)[:,:,end]
            energies[1:ns+1] = datatemp1[4]
            # append!(energies,datatemp1[4]) # Store the first batch of energies.
        
        
            # Store the output generated in the first simulation.
            pathnameaux = pathname*"/"*string(l)
            writedlm(pathnameaux*"_1_1.csv",pulledindicesT,',')
            writedlm(pathnameaux*"_1_2.csv",Int.(dirsT)) # Save only the numbers, not the whole enum type.
            writedlm(pathnameaux*"_1_3.csv",newcoordsT,',')
        
            # For each of the remaining temperatures, employ the Metropolis-Hastings algorithm to store the information about the
            # visited configurations at the current temperature.
            cont = ns+2
            for k in 2:length(temperatures)
                temp = reverse(temperatures)[k] # temperature at which the simulation is performed
                proteinaux = Protein(laststate,protein.HPlist,protein.geometry) # Protein structure for the simulation.
                pulledindicesT,dirsT,newcoordsT,enstates = metropolis(N,nums,temp,proteinaux,pfmodel)  # Perfom the simulation.
                # append!(energies,enstates) # Store the visited energies.
                energies[cont:cont+(ns-1)] = enstates[2:end] # Store the visited energies.
                cont = cont+ns 
                
                
                # Store the output generated in the first simulation.
                st = "_"*string(k)
                writedlm(pathnameaux*st*"_1.csv",pulledindicesT,',')
                writedlm(pathnameaux*st*"_2.csv",Int.(dirsT),',')
                writedlm(pathnameaux*st*"_3.csv",newcoordsT,',')
                
                laststate = reconstructStates(N,laststate,protein.HPlist,pulledindicesT,dirsT,newcoordsT,protein.geometry)[:,:,end] 
                println("Temperature Progress: $k /$nTs")
            end
            writedlm(pathname*"/energiesFolder/"*string(l)*"_energies.csv",energies,',') # Save all of the visted energies.
            println("Number of runs progress : $l /$nruns", )
        end
        
        
    elseif geometry == square2D
        # Now, perform a Metropolis-Hastings simulation for each temperature in `temperatures`. Sweep the tempeartures `nruns` times.
        @sync @distributed for l in 1:nruns
            energies = zeros(Float64,Int((ns*length(temperatures))+1)) # Energies for the current run.

            # Perform the first simulation.
            datatemp1 = metropolis(N,nums,temperatures[end],protein,pfmodel)
            pulledindicesT = datatemp1[1] # The next three variables will be used to perform subsequent simulations.
            dirsT = datatemp1[2]
            newcoordsT = datatemp1[3]
            laststate = reconstructStates(N,protein.edo,protein.HPlist,pulledindicesT,dirsT,newcoordsT,protein.geometry)[:,:,end]
            energies[1:ns+1] = datatemp1[4] # Store the first batch of energies.
    
            # Store the output generated in the first simulation.
            pathnameaux = pathname*"/"*string(l)
            writedlm(pathnameaux*"_1_1.csv",pulledindicesT,',')
            writedlm(pathnameaux*"_1_2.csv",Int.(dirsT),',') # Save only the Int value, not the whole enum info.
            writedlm(pathnameaux*"_1_3_1.csv",newcoordsT[1],',')
            writedlm(pathnameaux*"_1_3_2.csv",newcoordsT[2],',')

    
            # For each of the remaining temperatures, employ the Metropolis-Hastings algorithm to store the information about the
            # visited configurations at the current temperature.
            cont = ns+2
            for k in 2:length(temperatures)
                temp = reverse(temperatures)[k] # temperature at which the simulation is performed
                proteinaux = Protein(laststate,protein.HPlist,protein.geometry) # Auxiliary protein structure for the simulation.
                pulledindicesT,dirsT,newcoordsT,enstates = metropolis(N,nums,temp,proteinaux,pfmodel)  # Perform the simulation.
                energies[cont:cont+(ns-1)] = enstates[2:end] # Store the visited energies.
                cont = cont+ns 

                # Store the output generated in the first simulation.
                st = "_"*string(k)
                writedlm(pathnameaux*st*"_1.csv",pulledindicesT,',')
                writedlm(pathnameaux*st*"_2.csv",Int.(dirsT),',')
                writedlm(pathnameaux*st*"_3_1.csv",newcoordsT[1],',')
                writedlm(pathnameaux*st*"_3_2.csv",newcoordsT[2],',')
                
                laststate = reconstructStates(N,laststate,protein.HPlist,pulledindicesT,dirsT,newcoordsT,protein.geometry)[:,:,end] 
                println("Temperature Progress: $k /$nTs")
            end

            writedlm(pathname*"/energiesFolder/"*string(l)*"_energies.csv",energies,',') # Save all of the visted energies.
            println("Number of runs progress : $l /$nruns", )

        end

    end
    
    println("Simulation succesfully completed")
end
