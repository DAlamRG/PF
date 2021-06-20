using Plots
using Base: Tuple
# Here I do a couple of tests for my functions.
using DelimitedFiles

include("./HP-model2D.jl")
include("./visualizeHP2D.jl")
gr()


testprotein=Protein2D([[20 10];[20 11];[20 12];[20 13];[20 14];[20 15];[20 16];[20 17];[20 18];[20 19];[20 20];[20 21];[20 22];[20 23];
[20 24];[20 25];[20 26];[20 27];[20 28];[20 29]],[-1,1,-1,1,1,-1,-1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,-1],square2D)


#testproteinT=Protein2D([[7 2];[7 3];[7 4];[7 5];[7 6];[7 7];[7 8];[7 9];[7 10];[7 11]],[1,-1,-1,1,-1,1,-1,1,-1,-1],triangular2D)
# infoT=HP2Dmet(20,4000,0.2,testproteinT)
# reconT=reconstructStates2D(20,testproteinT.edo,testproteinT.HPlist,infoT[2],infoT[3],infoT[4],testproteinT.geometry)

infosquare2D=mainHP2Dmet(40,1000,0.1,3,21,testprotein) # Employ metropolis to simulate the protein folding process across a range of temperatures.
# Save the generated states.

writedlm("/Users/pedroruiz/Desktop/Diego/PF/Data/HPMet2D1_1.",infosquare2D[1])
writedlm("/Users/pedroruiz/Desktop/Diego/PF/Data/HPMet2D1_2.",infosquare2D[2])
writedlm("/Users/pedroruiz/Desktop/Diego/PF/Data/HPMet2D1_3.",infosquare2D[3])

visHP2D(infosquare2D[end],testprotein.HPlist,40,40,testprotein.geometry)
