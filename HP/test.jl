

# This script employs the main function to perform aand save the simulations.


include("./HP-model.jl")
# include("./HP_WL.jl")



# testProteinfcc = Protein([[7 2 9];[7 3 9];[7 4 9];[7 5 9];[7 6 9];[7 7 9];[7 8 9];[7 9 9];[7 10 9];[7 11 9];[7 12 9];[7 13 9];[7 14 9];[7 15 9];[7 16 9]],[H,P,N,h,H,h,P,X,P,h,N,h,H,P,H],fcc)

# testProtein=Protein([[10 5];[10 6];[10 7];[10 8];[10 9];[10 10];[10 11];[10 12];[10 13]],[P,H,N,H,P,X,H,H,P],triangular2D

seq_64 = Amin[H,H,H,H,H,H,H,H,H,H,H,H,P,H,P,H,P,P,H,H,P,P,H,H,P,P,H,P,P,H,H,P,P,H,H,P,P,H,P,P,H,H,P,P,H,H,P,P,H,P,H,P,H,H,H,H,H,H,H,H,H,H,H,H]

# seq_20 = Amin[]

testProtein64 = Protein(hcat(Int16[10 for i in 1:64],Int16(1+6):Int16(64+6)),seq_64,square2D)


# display(@time main_met_1(80,20,0.01,1.0,4,10,testProtein64,HP1_model,"simu8"))
# display(@time main_met(80,100,0.01,1.0,100,10,testProtein64,HP1_model,"simu3"))
# display(@time wang_landau(80,testProtein64,200,HP1_model,"simu2")) 




