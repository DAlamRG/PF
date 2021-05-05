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
    plot!(edo[:,2],edo[:,1],lw=2,color="purple",label="")
    display(plt)
end
