using StatsBase: minimum, maximum
using LinearAlgebra: convert, vcat

# This script employs the main function to perform and save simulations.



#include("./HP_MET.jl")
include("./HP_WL.jl")











# First, define a dictionary to translate between aminoacid symbols (from 1-letter to 3-letter). 
const dict1_3 =  Dict{String,Amin}("R" => ARG, "H" => HIS, "K" => LYS, "D" => ASP, "E" => GLU,
"N" => ASN, "C" => CYS, "Q" => GLN, "S" => SER, "T" => THR, "Y" => TYR, "A" => ALA, "G" => GLY,
"I" => ILE, "L" => LEU, "M" => MET, "F" => PHE, "P" => PRO, "W" => TRP, "V" => VAL)

"""
    convert_Amin(str)
Given a string of characters encoding the protein's amino acid sequence `str`; returns the equivalent sequence as an `Amin` array.
"""
function convert_Amin(str::String)
    split_str = split(str,"") # Split the string into characters.
    HPlist = fill(CYS,length(split_str)) # Declare a vector of type AMIN.
    for i in 1:length(split_str)
       HPlist[i] = dict1_3[split_str[i]] 
    end
    return HPlist
end

# Now, we define some sequences.
# A good website to search for sequences is https://www.rcsb.org/ or https://www.uniprot.org

seq_64 = Amin[H,H,H,H,H,H,H,H,H,H,H,H,P,H,P,H,P,P,H,H,P,P,H,H,P,P,H,P,P,H,H,P,P,H,H,P,P,H,P,P,H,H,P,P,H,H,P,P,H,P,H,P,H,H,H,H,H,H,H,H,H,H,H,H]

hemo = convert_Amin("MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFFDSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTFAQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQASWQKMVTGVASALSSRYH")

ribonuclease = convert_Amin("MALEKSLVRLLLLVLILLVLGWVQPSLGKESRAKKFQRQHMDSDSSPSSSSTYCNQMMRRRNMTQGRCKPVNTFVHEPLVDVQNVCFQEKVTCKNGQGNCYKSNSSMHITDCRLTNGSRYPNCAYRTSPKERHIIVACEGSPYVPVHFDASVEDST")

crambin = convert_Amin("TTCCPSIVARSNFNVCRLPGTPEALCATYTGCIIIPGATCPGDYAN")

chignolin = convert_Amin("YYDPETGTWY") # https://www.rcsb.org/structure/5AWL

trp_cage = convert_Amin("DAYAQWLKDGGPSSGRPPPS") # https://www.rcsb.org/structure/2JOF

villin = convert_Amin("LSDEDFKAVFGMTRSAFANLPLWLQQHLLKEKGLF") # https://www.rcsb.org/structure/2F4K # Article includes NLE aminoacid!

α3D = convert_Amin("MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAYKGKGNPEVEALRKEAAAIRDELQAYRHN") # https://www.rcsb.org/structure/2A3D













chignolin_Full1_square = Protein(hcat(Int16[10 for i in 1:10],Int16(1+6):Int16(10+6)), translate_HPlist(chignolin,Full1,true),square2D)


#trp_cage_HP1_triangular = Protein(hcat(Int16[10 for i in 1:20],Int16(1+6):Int16(20+6)), translate_HPlist(trp_cage,HP1,true),triangular2D)
#trp_cage_HPNX_triangular = Protein(hcat(Int16[10 for i in 1:20],Int16(1+6):Int16(20+6)), translate_HPlist(trp_cage,HPNX,true),triangular2D)
#trp_cage_Full1_triangular = Protein(hcat(Int16[10 for i in 1:20],Int16(1+6):Int16(20+6)), translate_HPlist(trp_cage,Full1,true),triangular2D)

#trp_cage_HP1_fcc = Protein(hcat(Int16[10 for i in 1:20],Int16(1+6):Int16(20+6),Int16[10 for i in 1:20]),translate_HPlist(trp_cage,HP1,true),fcc)
#trp_cage_HPNX_fcc = Protein(hcat(Int16[10 for i in 1:20],Int16(1+6):Int16(20+6),Int16[10 for i in 1:20]),translate_HPlist(trp_cage,HPNX,true),fcc)
#trp_cage_Full1_fcc = Protein(hcat(Int16[10 for i in 1:20],Int16(1+6):Int16(20+6),Int16[10 for i in 1:20]),translate_HPlist(trp_cage,Full1,true),fcc)







# "name" should follow the format "WL_chignolin_HP1_square"

#display(@time main_met(22,1000,0.01,1.0,600,16,chignolin_HP1_square,HP1_model,"MET_chignolin_HP1_square"))

display(@time wang_landau(34,trp_cage_Full1_triangular,85,4,6,Full1_model,"WL_trp_cage_Full1_triangular"))













