# I write a function which returns a plot which allows for better visualization of the amino acids type, positions and bonds.
using Plots




# Write the visaulization function. 
"""
    visHP2D(red,edo)

Given the array `red` containing the protein´s configuration, a matrix encoding the aminoacids positions `edo`; returns 
a heat map with the protein´s configuration.
"""
function visHP2D(red,edo)
    N=size(red)[1]
    plt=heatmap(1:N,1:N,red,xlabel="",ylabel="",label="",title="Current protein configuration")
    plot!(edo[:,2],edo[:,1],lw=2,markershape=:circle,markersize=3,color="purple",label="")
    display(plt)
end











# I would also like to be able to visulaize the protein´s energies across iterations of the Metropolis scheme.
"""
    energiesHP2D(T,energies)

Given the temperature `T`, and an array containing the visited energies during the simulation `energies`; returns 
a plot of the energies as a function of iteration step.
"""
function energiesHP2D(T,energies)
    pp=plot(1:length(energies),energies,color="green",lw=2,alpha=0.7,xlabel="iteration step",
    ylabel="E(T)",label="",title="Visited energies at T=$T")
    display(pp)
end




