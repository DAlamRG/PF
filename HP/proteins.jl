# This script contains some protein chains. 


include("./HP-pullmoves.jl")

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

chignolin = convert_Amin("YYDPETGTWY") # https://www.rcsb.org/structure/5AWL

trp_cage = convert_Amin("DAYAQWLKDGGPSSGRPPPS") # https://www.rcsb.org/structure/2JOF

villin = convert_Amin("LSDEDFKAVFGMTRSAFANLPLWKQQHLKKEKGLF") # https://www.rcsb.org/fasta/entry/1YRF/display # Article includes NLE aminoacid!

α3D = convert_Amin("MGSWAEFKQRLAAIKTRLQALGGSEAELAAFEKEIAAFESELQAYKGKGNPEVEALRKEAAAIRDELQAYRHN") # https://www.rcsb.org/structure/2A3D