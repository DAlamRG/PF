using Plots: length
using Plots
using Base: Tuple
# Here I do a couple of tests for my functions.
using DelimitedFiles

include("./HP-model2D.jl")
include("./visualizeHP2D.jl")
gr()


testprotein=Protein2D([[20 10];[20 11];[20 12];[20 13];[20 14];[20 15];[20 16];[20 17];[20 18];[20 19];[20 20];[20 21];[20 22];[20 23];
[20 24];[20 25];[20 26];[20 27];[20 28];[20 29]],[-1,1,-1,1,1,-1,-1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,-1],square2D)


# testproteinTriang=Protein2D([[7 7];[7 8];[7 9];[7 10];[7 11];[7 12];[7 13];[7 14];[7 15];[7 16];[7 17];[7 18]],[1,-1,-1,1,-1,1,-1,1,-1,-1,1,-1],triangular2D)
# infoT=HP2Dmet(20,4000,0.2,testproteinT)
# reconT=reconstructStates2D(20,testproteinT.edo,testproteinT.HPlist,infoT[2],infoT[3],infoT[4],testproteinT.geometry)



 mainHP2Dmet(38,5,0.1,1.1,10,2,testprotein,"trial2DSQUARE1")



# Now that we have a whole simulation over a range of temperatures, I make a little animation.""