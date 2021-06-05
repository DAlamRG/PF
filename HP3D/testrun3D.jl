include("HP-model3D.jl")

plotly()



testProteinfcc=Protein([[7 2 9];[7 3 9];[7 4 9];[7 5 9];[7 6 9];[7 7 9];[7 8 9];[7 9 9];[7 10 9];[7 11 9];[7 12 9];[7 13 9];[7 14 9]],
[1,-1,-1,1,-1,1,-1,1,-1,-1,1,-1,-1],fcc)

infofcc=HP3Dmet(20,1000,0.3,testProteinfcc) # Employ metropolis to simulate the protein folding process.

recon=reconstructStates3D(20,testProteinfcc.edo,testProteinfcc.HPlist,infofcc[2],infofcc[3],infofcc[4],testProteinfcc.geometry) # Reconstruct the visited states during the simulation.

display(infofcc[1])

visHP3D(recon[:,:,end],testProteinfcc.HPlist,20,testProteinfcc.geometry)