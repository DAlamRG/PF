using StatsBase: minimum, maximum
using LinearAlgebra: convert, vcat

# This script employs the main function to perform and save simulations.




# include("./HP-model.jl")
include("./HP_WL.jl")


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



# First, define a `Protein` structure:
# A good website to search for sequences is uniprot.org

# seq_64 = Amin[H,H,H,H,H,H,H,H,H,H,H,H,P,H,P,H,P,P,H,H,P,P,H,H,P,P,H,P,P,H,H,P,P,H,H,P,P,H,P,P,H,H,P,P,H,H,P,P,H,P,H,P,H,H,H,H,H,H,H,H,H,H,H,H]
# testProtein64 = Protein(hcat(Int16[10 for i in 1:64],Int16(1+6):Int16(64+6)),seq_64,square2D)

# seq_hemo = convert_Amin("MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQASWQKMVTGVASALSSRYH")
# seq_hemo = translate_HPlist(seq_hemo,HP1,true)

# seq_ribonuclease = convert_Amin("MALEKSLVRLLLLVLILLVLGWVQPSLGKESRAKKFQRQHMDSDSSPSSSSTYCNQMMRRRNMTQGRCKPVNTFVHEPLVDVQNVCFQEKVTCKNGQGNCYKSNSSMHITDCRLTNGSRYPNCAYRTSPKERHIIVACEGSPYVPVHFDASVEDST")

seq_crambinaux = convert_Amin("TTCCPSIVARSNFNVCRLPGTPEALCATYTGCIIIPGATCPGDYAN")
seq_crambinaux = translate_HPlist(seq_crambinaux,Full1,true)
crambinaux_462D = Protein(hcat(Int16[20 for i in 1:46],Int16(1+6):Int16(46+6)),seq_crambinaux,triangular2D)

# seq_crambin = Amin[P,P,H,H,P,P,H,H,H,P,P,P,H,P,H,H,P,H,P,P,P,P,P,H,H,H,H,P,H,P,P,H,H,H,H,P,P,H,P,H,P,P,P,H,H,P]

# crambin_462D = Protein(hcat(Int16[20 for i in 1:46],Int16(1+6):Int16(46+6)),seq_crambin,triangular2D)
# crambin_463D = Protein(hcat(Int16[20 for i in 1:46],Int16(1+6):Int16(46+6),Int16[30 for i in 1:46]),seq_crambin,fcc)


# display(@time main_met(62,130,0.01,1.2,85,15,crambin_462D,HP1_model,"simu2"))
display(@time wang_landau(62,crambinaux_462D,210,Full1_model,"simu6"))


# display(@time main_met(62,112,0.01,6.0,72,16,crambinaux_462D,HP2_model,"simu11"))
# display(@time main_met(62,115,0.01,8.0,75,16,crambinaux_462D,hHPNX_model,"simu12"))
# display(@time main_met(62,120,0.01,6.0,70,16,crambinaux_462D,YhHX_model,"simu13"))
# display(@time main_met(62,125,0.01,6.0,70,16,crambinaux_462D,Full1_model,"simu14"))











