using Base: Int8
using DelimitedFiles
include("HP-model3D.jl")

plotly()



testProteinfcc=Protein([[7 2 9];[7 3 9];[7 4 9];[7 5 9];[7 6 9];[7 7 9];[7 8 9];[7 9 9];[7 10 9];[7 11 9];[7 12 9];[7 13 9];[7 14 9]],
[1,-1,-1,1,-1,1,-1,1,-1,-1,1,-1,-1],fcc)

# Employ metropolis to simulate the protein folding process across a range of temperatures.

 mainHP3Dmet(20,5,0.1,3,11,3,testProteinfcc,"test2") 

 