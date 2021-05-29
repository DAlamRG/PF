

using Plots
    
    
    
    
# Write the visualization function. 
"""
    visHP3D(edo,HPlist,xlims,ylims,zlims)
Given a matrix encoding the aminoacids positions `edo`, and the aminoacid sequence `HPlist`; returns 
a plot displaying the protein configuration.
"""
function visHP3D(edo,HPlist,xlims,ylims,zlims)
    colours=[]
    for i in 1:length(HPlist)
        el=HPlist[i]
        if el == 1
            push!(colours,"green")
        else
            push!(colours,"red")
        end
    end

    
    plt=plot(edo[:,2],edo[:,1],edo[:,3],lw=2,markershape=:circle,markercolor=colours ,markersize=7,color="purple",label="",xlabel="",
    ylabel="",title="Current protein configuration ")
    xlims!(xlims)
    ylims!(ylims)
    zlims!(zlims)
    display(plt)
end
    
    
    
    
    
    
    
    
    
    
    
    
# I would also like to be able to visulaize the protein´s energies across iterations of the Metropolis scheme.
    
"""
        energiesHP(T,energies)
    
Given the temperature `T`, and an array containing the visited energies during the simulation `energies`; returns 
a plot of the energies as a function of iteration step.
"""
function energiesHP(T,energies)
    pp=plot(1:length(energies),energies,color="green",lw=2,alpha=0.7,xlabel="iteration step",
    ylabel="E",label="",title="Visited energies at T=$T")
    display(pp)
end