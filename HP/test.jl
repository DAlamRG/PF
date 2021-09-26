using LinearAlgebra: convert

# This script employs the main function to perform and save the simulations.
include("./HP-model.jl")





# First, define a `Protein` structure:
# A good website to search for sequences is uniprot.org


# testProteinfcc = Protein([[7 2 9];[7 3 9];[7 4 9];[7 5 9];[7 6 9];[7 7 9];[7 8 9];[7 9 9];[7 10 9];[7 11 9];[7 12 9];[7 13 9];[7 14 9];[7 15 9];[7 16 9]],[H,P,N,h,H,h,P,X,P,h,N,h,H,P,H],fcc)
# testProtein = Protein([[10 5];[10 6];[10 7];[10 8];[10 9];[10 10];[10 11];[10 12];[10 13]],[P,H,N,H,P,X,H,H,P],triangular2D

# seq_64 = Amin[H,H,H,H,H,H,H,H,H,H,H,H,P,H,P,H,P,P,H,H,P,P,H,H,P,P,H,P,P,H,H,P,P,H,H,P,P,H,P,P,H,H,P,P,H,H,P,P,H,P,H,P,H,H,H,H,H,H,H,H,H,H,H,H]
# testProtein64 = Protein(hcat(Int16[10 for i in 1:64],Int16(1+6):Int16(64+6)),seq_64,square2D)

# Here I define the protein corresponding to  Hemoglobin subunit gamma-2
const dict1_3 =  Dict{String,Amin}("R" => ARG, "H" => HIS, "K" => LYS, "D" => ASP, "E" => GLU,
"N" => ASN, "C" => CYS, "Q" => GLN, "S" => SER, "T" => THR, "Y" => TYR, "A" => ALA, "G" => GLY,
"I" => ILE, "L" => LEU, "M" => MET, "F" => PHE, "P" => PRO, "W" => TRP, "V" => VAL)


"""
    convert_Amin(str)
Given a string of characters encoding the protein's amino acid sequence `str`; returns the equivalent sequence as an `Amin` array.
"""
function convert_Amin(str::String)
    split_str = split(str,"")
    HPlist = fill(CYS,length(split_str))
    for i in 1:length(split_str)
       HPlist[i] = dict1_3[split_str[i]] 
    end
    return HPlist
end

seq_hemo = convert_Amin("MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQASWQKMVTGVASALSSRYH")
seq_hemo = translate_HPlist(seq_hemo,HP1,true)

# seq_ribonuclease = convert_Amin("MALEKSLVRLLLLVLILLVLGWVQPSLGKESRAKKFQRQHMDSDSSPSSSSTYCNQMMRRRNMTQGRCKPVNTFVHEPLVDVQNVCFQEKVTCKNGQGNCYKSNSSMHITDCRLTNGSRYPNCAYRTSPKERHIIVACEGSPYVPVHFDASVEDST")
seq_crambin = convert_Amin("TTCCPSIVARSNFNVCRLPGTPEALCATYTGCIIIPGATCPGDYAN")
seq_crambin = translate_HPlist(seq_crambin,HP1,true)


crmabin_46 = Protein(hcat(Int16[20 for i in 1:46],Int16(1+6):Int16(46+6)),seq_crambin,fcc)


display(@time main_met_1(160,80,0.01,1.0,50,1,hemoglobin147,HP1_model,"simu6"))

# display(@time main_met(38,100,0.01,1.0,30,1,hemoglobin21,HP1_model,"simu6"))
# display(@time main_met(55,100,0.01,1.0,30,1,hemoglobin42,HP1_model,"simu7"))
# display(@time main_met(75,100,0.01,1.0,30,1,hemoglobin63,HP1_model,"simu8"))
# display(@time main_met(96,100,0.01,1.0,30,1,hemoglobin84,HP1_model,"simu9"))
# display(@time main_met(120,100,0.01,1.0,30,1,hemoglobin105,HP1_model,"simu10"))
# display(@time main_met(140,100,0.01,1.0,30,1,hemoglobin126,HP1_model,"simu11"))
# display(@time main_met(160,100,0.01,1.0,30,1,hemoglobin147,HP1_model,"simu12"))





