# Here I do a couple of tests for my functions.

include("./HP-model2D.jl")
gr()


testprotein=Protein2D([[20 10];[20 11];[20 12];[20 13];[20 14];[20 15];[20 16];[20 17];[20 18];[20 19];[20 20];[20 21];[20 22];[20 23];
[20 24];[20 25];[20 26];[20 27];[20 28];[20 29]],[-1,1,-1,1,1,-1,-1,1,-1,1,1,-1,1,-1,-1,1,1,-1,1,-1],square2D)

info=HP2Dmet(40,4000,0.2,testprotein)
recon=reconstructStates2D(40,testprotein.edo,testprotein.HPlist,info[2],info[3],info[4],testprotein.geometry)




testproteinT=Protein2D([[7 2];[7 3];[7 4];[7 5];[7 6];[7 7];[7 8];[7 9];[7 10];[7 11]],[1,-1,-1,1,-1,1,-1,1,-1,-1],triangular2D)

# infoT=HP2Dmet(20,4000,0.2,testproteinT)
# reconT=reconstructStates2D(20,testproteinT.edo,testproteinT.HPlist,infoT[2],infoT[3],infoT[4],testproteinT.geometry)




println("Reconstruction of square was= ",info[end] == recon )
# println("Reconstruction of triangle was= ",infoT[end] == reconT )