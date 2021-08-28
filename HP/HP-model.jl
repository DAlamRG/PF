using Base: Int16
# This script runs a Metropolis simulation, and stores the data.

using DelimitedFiles

include("./Energy.jl")











# Now that I am able to generate a new state by performing a pull move, and knowing how to compute the energy of any
# given configuration, I proceed to implement a simulation using the Metropolis-Hastings algorithm.

"""
    metropolis(N,nums,T,protein,pfmodel) 

Given a 2D/3D array size `N`, a number of Monte-Carlo sweeps `nums`, the temperature T, a protein structure `protein` encoding the protein´s 
configurationand a aminoacid interaction model `pfmodel`; returns the final configuration and the visited energies after performing a simulation 
using the Metropolis-Hastings algorithm.
"""
function metropolis(N,nums,T,protein,pfmodel)  

    β=1/T
    edo=protein.edo
    HPlist=protein.HPlist
    geometry=protein.geometry
    ns=nums*length(HPlist) # Total number of iterations


    states=ones(Int16,(length(HPlist),length(edo[1,:]),ns+1)) # `states` stores the amino acids´ coordinates for each visited configuration.
    states[:,:,1]=edo # First coordinates are those of the input state.
    
    difes=zeros(Float64,ns) # Stores the energy difference ΔH between states.
    enstates=zeros(Float64,ns+1) # Stores the energy of each state visited during the simulation.
    enstates[1]=energy(N,edo,HPlist,geometry,pfmodel) # Compute the energy of the initial state and store it.
    
    npullstates=zeros(Float64,ns+1) # Stores the number of pull moves for each configuration.
    npullstates[1]=countpull(N,edo,HPlist,geometry)[1] # Count the number of pull moves for the initial state and store it.

        
    # Now, apply Metropolis-Hastings.
    
    if geometry == fcc || geometry == triangular2D
        pulledindices=zeros(Int16,ns) # Stores the moved indices. If no index is moved, the value is zero.
        newcoords=zeros(Int16,(ns,length(edo[1,:]))) # Stores the new coordinate for the pulled index. If the configuration stays the same, the new coordinate is just the original position.
        dirs=directions[] # Store the direction in which the chain was pulled. If no pull move was peroformed, the direction is `nonetaken`.

        for l in 2:ns+1 
        
            # Generate a new state.
            newedo,totalpull,indpulled,newcoord,dirpulled = pullMove(N,states[:,:,l-1],HPlist,geometry)
            newenergy=energy(N,newedo,HPlist,geometry,pfmodel) # Compute the energy after the pull move.
            
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
                    newcoords[l-1,:]=repeat([0],length(edo[1,:]))
                    push!(dirs,nonetaken)
                end
            end
        end

        return (pulledindices,dirs,newcoords,enstates)
        

    elseif geometry == square2D
        pulledindices=zeros(Int8,ns) # Stores the moved indices. If no index is moved, the value is zero.
        newcoords1=zeros(Int8,(ns,length(edo[1,:]))) # Stores the new coordinate for the pulled index. 
        newcoords2=zeros(Int8,(ns,length(edo[1,:]))) # Stores the new coordinate for the momomer directly before/after the one being moved. 
        dirs=directions[] # Store the direction in which the chain was pulled. If no pull move was performed, the direction is `nonetaken`.

        #Generate new states by performing pull moves.
        for l in 2:ns+1 
            # Generate a new state.
            newedo,totalpull,indpulled,newcoord1,newcoord2,dirpulled = pullMove(N,states[:,:,l-1],HPlist,geometry)
            newenergy=energy(N,newedo,HPlist,geometry,pfmodel) # Compute the energy after the pull move.
            
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
                    newcoords1[l-1,:]=repeat([0],length(edo[1,:]))
                    newcoords2[l-1,:]=repeat([0],length(edo[1,:]))
                    push!(dirs,nonetaken)
                end
            end
        end

        newcoords=(newcoords1,newcoords2)
        return (pulledindices,dirs,newcoords,enstates)

    end

end





















 # Next, I write a function which performs multiple simulations over an array of temperatures, and saves the information generated.
 """
     main_met(N,nums,ti,tf,nTs,nruns,protein,pfmodel,name)

Given a 2D/3D array size `N`, a number of Monte-Carlo sweeps per temperature `nums`, an initial(final) temperature `ti(tf)`, the number of 
temperatures to be visited `nTs`, a number of independent temperature sweeps `nruns`, a structure `protein` encoding the protein´s 
configuration, an aminoacid interaction model `pfmodel`, and a name for the simulation output `name`; writes the information needed to 
recreate the visited states after performing a simulation using the Metropolis-Hastings algorithm.
"""
function main_met(N,nums,ti,tf,nTs,nruns,protein,pfmodel,name)
    temperatures=range(ti,stop=tf,length=nTs) # Decalre a range of temperatures to be visited.

    # Create the directory which will contain the dat collected trough the simulation.
    pathstring="./output/"
    pathname=pathstring*name
    mkdir(pathname)
    writedlm(pathname*"/temperatures.csv",temperatures,',')
    writedlm(pathname*"/initialconf.csv",protein.edo,',')
    writedlm(pathname*"/HPlist.csv",protein.HPlist,',')
    writedlm(pathname*"/latticesize.csv",N,',')
    writedlm(pathname*"/mc_sweeps.csv",nums,',')
    writedlm(pathname*"/geometry.csv",Int(protein.geometry),',')
    writedlm(pathname*"/pfmodel.csv",Int(pfmodel),',')

    geometry=protein.geometry

    if geometry == fcc || geometry == triangular2D
        # Now, perform a Metropolis-Hastings simulation for each temperature in `temperatures`. Sweep the tempeartures `nruns` times.
        for l in 1:nruns
            energies=Float64[] # Energies for the current run.
        
            # Perform the first simulation.
            datatemp1=metropolis(N,nums,temperatures[end],protein,pfmodel)
            pulledindicesT=datatemp1[1] # The next three variables will be used to perform subsequent simulations.
            dirsT=datatemp1[2]
            newcoordsT=datatemp1[3]
            laststate=reconstructStates(N,protein.edo,protein.HPlist,pulledindicesT,dirsT,newcoordsT,protein.geometry)[:,:,end]
            append!(energies,datatemp1[4]) # Store the first batch of energies.
        
        
            # Store the output generated in the first simulation.
            pathnameaux=pathname*"/"*string(l)
            writedlm(pathnameaux*"_1_1.csv",pulledindicesT,',')
            writedlm(pathnameaux*"_1_2.csv",Int.(dirsT)) # Save only the numbers, not the whole enum type.
            writedlm(pathnameaux*"_1_3.csv",newcoordsT,',')
        
            # For each of the remaining temperatures, employ the Metropolis-Hastings algorithm to store the information about the
            # visited configurations at the current temperature.
            for k in 2:length(temperatures)
                temp=reverse(temperatures)[k] # temperature at which the simulation is performed
                proteinaux=Protein(laststate,protein.HPlist,protein.geometry) # Protein structure for the simulation.
                pulledindices,dirs,newcoords,enstates=metropolis(N,nums,temp,proteinaux,pfmodel)  # Perfom the simulation.
                append!(energies,enstates) # Store the visited energies.
                
                pulledindicesT=pulledindices # The next three variables will be used to perform subsequent simulations.
                dirsT=dirs
                newcoordsT=newcoords
                
                # Store the output generated in the first simulation.
                st="_"*string(k)
                writedlm(pathnameaux*st*"_1.csv",pulledindicesT,',')
                writedlm(pathnameaux*st*"_2.csv",Int.(dirsT),',')
                writedlm(pathnameaux*st*"_3.csv",newcoordsT,',')
                
                laststate=reconstructStates(N,laststate,protein.HPlist,pulledindicesT,dirsT,newcoordsT,protein.geometry)[:,:,end] 
            end
            writedlm(pathnameaux*"_energies.csv",energies,',') # Save all of the visted energies.
        end
        
        
    elseif geometry == square2D
        # Now, perform a Metropolis-Hastings simulation for each temperature in `temperatures`. Sweep the tempeartures `nruns` times.
        for l in 1:nruns
            energies=Float64[] # Energies for the current run.

            # Perform the first simulation.
            datatemp1=metropolis(N,nums,temperatures[end],protein,pfmodel)
            pulledindicesT=datatemp1[1] # The next three variables will be used to perform subsequent simulations.
            dirsT=datatemp1[2]
            newcoordsT=datatemp1[3]
            laststate=reconstructStates(N,protein.edo,protein.HPlist,pulledindicesT,dirsT,newcoordsT,protein.geometry)[:,:,end]
            append!(energies,datatemp1[4]) # Store the first batch of energies.
    
    
            # Store the output generated in the first simulation.
            pathnameaux=pathname*"/"*string(l)
            writedlm(pathnameaux*"_1_1.csv",pulledindicesT,',')
            writedlm(pathnameaux*"_1_2.csv",Int.(dirsT),',') # Save only the Int value, not the whole enum info.
            writedlm(pathnameaux*"_1_3_1.csv",newcoordsT[1],',')
            writedlm(pathnameaux*"_1_3_2.csv",newcoordsT[2],',')

    
            # For each of the remaining temperatures, employ the Metropolis-Hastings algorithm to store the information about the
            # visited configurations at the current temperature.
            for k in 2:length(temperatures)
                temp=reverse(temperatures)[k] # temperature at which the simulation is performed
                proteinaux=Protein(laststate,protein.HPlist,protein.geometry) # Auxiliary protein structure for the simulation.
                pulledindices,dirs,newcoords,enstates=metropolis(N,nums,temp,proteinaux,pfmodel)  # Perform the simulation.
                append!(energies,enstates) # Store the visited energies.

                pulledindicesT=pulledindices # The next three variables will be used to perform subsequent simulations.
                dirsT=dirs
                newcoordsT=newcoords
                
                # Store the output generated in the first simulation.
                st="_"*string(k)
                writedlm(pathnameaux*st*"_1.csv",pulledindicesT,',')
                writedlm(pathnameaux*st*"_2.csv",Int.(dirsT),',')
                writedlm(pathnameaux*st*"_3_1.csv",newcoordsT[1],',')
                writedlm(pathnameaux*st*"_3_2.csv",newcoordsT[2],',')
                
                laststate=reconstructStates(N,laststate,protein.HPlist,pulledindicesT,dirsT,newcoordsT,protein.geometry)[:,:,end] 
            end


            writedlm(pathnameaux*"_energies.csv",energies,',') # Save all of the visted energies.

        end

    end
    
    println("Simulation succesfully completed")
end


