using Base: Int8
using DelimitedFiles
include("HP-model3D.jl")

plotly()



testProteinfcc=Protein([[7 2 9];[7 3 9];[7 4 9];[7 5 9];[7 6 9];[7 7 9];[7 8 9];[7 9 9];[7 10 9];[7 11 9];[7 12 9];[7 13 9];[7 14 9]],
[1,-1,-1,1,-1,1,-1,1,-1,-1,1,-1,-1],fcc)

infofcc=mainHP3Dmet(20,100,0.1,3,21,testProteinfcc) # Employ metropolis to simulate the protein folding process across a range of temperatures.
# Save the generated states.

writedlm("/Users/pedroruiz/Desktop/Diego/PF/Data/HPMet3D1_1.",infofcc[1])
writedlm("/Users/pedroruiz/Desktop/Diego/PF/Data/HPMet3D1_2.",infofcc[2])
writedlm("/Users/pedroruiz/Desktop/Diego/PF/Data/HPMet3D1_3.",infofcc[3])


visHP3D(infofcc[end],testProteinfcc.HPlist,20,testProteinfcc.geometry)
