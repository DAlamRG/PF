

# This script employs the main function to perform aand save the simulations.


include("./HP-model.jl")



testProteinfcc = Protein([[7 2 9];[7 3 9];[7 4 9];[7 5 9];[7 6 9];[7 7 9];[7 8 9];[7 9 9];[7 10 9];[7 11 9];[7 12 9];[7 13 9];[7 14 9];[7 15 9];[7 16 9]],[H,P,N,h,H,h,P,X,P,h,N,h,H,P,H],fcc)

# testProtein=Protein([[10 5];[10 6];[10 7];[10 8];[10 9];[10 10];[10 11];[10 12];[10 13]],[P,H,N,H,P,X,H,H,P],triangular2D)


main_met(19,10,0.01,0.6,12,5,testProteinfcc,hHPNX_model,"simu1")




