# Generalized pull-moves for both dimensions and the three geometries.

# First, declare the relevant structures for both dimensions.

using LinearAlgebra



@enum geometries::Int8 begin # First, define the type of geometries as an enum type.
    square2D = 1
    triangular2D = 2
    cubic = 3
    fcc = 4
end




"""
    geom_int(val)
Given an `Int` type value for the geometry; returns the equivalent `enum` value.
"""
function geom_int(val::Int)
    enumval = square2D
    if val == 1
        enumval = square2D
    elseif val == 2
        enumval = triangular2D
    elseif val == 3
        enumval = cubic
    elseif val == 4
        enumval = fcc
    end
    return enumval
end














# Declare the type of aminoacids available for this script.
@enum Amin::Int8 begin
    h = 1
    H = 2
    P = 3
    N = 4
    X = 5
    Y = 6
    CYS = 7
    MET = 8
    PHE = 9
    ILE = 10
    LEU = 11
    VAL = 12
    TRP = 13
    TYR = 14 
    ALA = 15
    GLY = 16
    THR = 17
    SER = 18
    ASN = 19
    GLN = 20
    ASP = 21
    GLU = 22 
    HIS = 23 
    ARG = 24
    LYS = 25
    PRO = 26
end

# Dictionary to translate between Int and aminoacid name.
const amin_dict = Dict{Int8,Amin}(1 => h, 2 => H, 3 => P, 4 => N, 5 => X, 6 => Y, 
7 => CYS, 8 => MET, 9 => PHE, 10 => ILE, 11 => LEU, 12 => VAL, 13 => TRP,
14 => TYR, 15 => ALA, 16 => GLY, 17 => THR, 18 => SER, 19 => ASN, 20 => GLN,
21 => ASP, 22 => GLU, 23 => HIS, 24 => ARG, 25 => LYS, 26 => PRO)















# Declare the different interaction model names.
@enum PFmodelname::Int8 begin
    HP1 = 1
    HP2 = 2
    HP3 = 3
    HPNX = 4
    hHPNX = 5
    YhHX = 6
    Full1 = 7
    Full2 = 8
end 

# Dictionary to translate between Int and model name.
const pfname_dict = Dict{Int8,PFmodelname}(1 => HP1, 2 => HP2, 3 => HP3, 4 => HPNX, 5 => hHPNX, 6 => YhHX, 7 => Full1, 8 => Full2)













# Dictionaries to assign positions.
const DictHP = Dict{Amin,Int8}(H => 1, P => 2, N => 3, X => 4, h => 5)
const DictYhHX = Dict{Amin,Int8}(Y => 1, h => 2, H => 3, X => 4)
const DictFull1 = Dict{Amin,Int8}(CYS => 1, MET => 2, PHE => 3, ILE => 4, LEU => 5,
VAL => 6, TRP => 7, TYR => 8, ALA => 9, GLY => 10, THR => 11, SER => 12, ASN => 13,
GLN => 14, ASP => 15, GLU => 16, HIS => 17, ARG => 18, LYS => 19, PRO => 20)

const HPNXintMatrix = [
    [-4 0 0 0] ; 
    [0 1 -1 0] ; 
    [0 -1 1 0] ; 
    [0 0 0 0]
]

const HP1intMatrix = [
    [-1 0] ; 
    [0 0] 
]

const HP2intMatrix = [
    [-3 -1] ; 
    [-1 0] 
]

const HP3intMatrix = [
    [-2.3 -1] ; 
    [-1 0]  
]

const hHPNXintMatrix = [
    [-3 0 0 0 -4] ; 
    [0 1 -1 0 0] ; 
    [0 -1 1 0 0] ; 
    [0 0 0 0 0] ; 
    [-4 0 0 0 2] 
]

const YhHXintMatrix = [
    [0 -1 -1 2] ; 
    [-1 2 -4 2] ; 
    [-1 -4 -3 0] ; 
    [2 2 0 0] 
]

const Full1intMatrix = [
    [-3.477 -2.240 -2.424 -2.410 -2.343 -2.258 -2.080 -1.892 -1.700 -1.101 -1.243 -1.306 -0.788 -0.835 -0.616 -0.179 -1.499 -0.771 -0.112 -1.196] ; 
    [-2.240 -1.901 -2.304 -2.286 -2.208 -2.079 -2.090 -1.834 -1.517 -0.897 -0.999 -0.893 -0.658 -0.720 -0.409 -0.209 -1.252 -0.611 -0.146 -0.788] ; 
    [-2.424 -2.304 -2.467 -2.530 -2.491 -2.391 -2.286 -1.963 -1.750 -1.034 -1.237 -1.178 -0.790 -0.807 -0.482 -0.419 -1.330 -0.805 -0.270 -1.076] ; 
    [-2.410 -2.286 -2.530 -2.691 -2.647 -2.568 -2.303 -1.998 -1.872 -0.885 -1.360 -1.037 -0.669 -0.778 -0.402 -0.439 -1.234 -0.854 -0.253 -0.991] ;
    [-2.343 -2.208 -2.491 -2.647 -2.501 -2.447 -2.222 -1.919 -1.728 -0.767 -1.202 -0.959 -0.524 -0.729 -0.291 -0.366 -1.176 -0.758 -0.222 -0.771] ;
    [-2.258 -2.079 -2.391 -2.568 -2.447 -2.385 -2.097 -1.790 -1.731 -0.756 -1.240 -0.933 -0.673 -0.642 -0.298 -0.335 -1.118 -0.664 -0.200 -0.886] ;
    [-2.080 -2.090 -2.286 -2.303 -2.222 -2.097 -1.867 -1.834 -1.565 -1.142 -1.077 -1.145 -0.884 -0.997 -0.613 -0.624 -1.383 -0.912 -0.391 -1.278] ;
    [-1.892 -1.834 -1.963 -1.998 -1.919 -1.790 -1.834 -1.335 -1.318 -0.818 -0.892 -0.859 -0.670 -0.687 -0.631 -0.453 -1.222 -0.745 -0.349 -1.067] ;
    [-1.700 -1.517 -1.750 -1.872 -1.728 -1.731 -1.565 -1.318 -1.119 -0.290 -0.717 -0.607 -0.371 -0.323 -0.235 -0.039 -0.646 -0.327 0.196 -0.374] ;
    [-1.101 -0.897 -1.034 -0.885 -0.767 -0.756 -1.142 -0.818 -0.290 0.219 -0.311 -0.261 -0.230 0.033 -0.097 0.443 -0.325 -0.050 0.589 -0.042] ;
    [-1.243 -0.999 -1.237 -1.360 -1.202 -1.240 -1.077 -0.892 -0.717 -0.311 -0.617 -0.548 -0.463 -0.342 -0.382 -0.192 -0.720 -0.247 0.155 -0.222] ;
    [-1.306 -0.893 -1.178 -1.037 -0.959 -0.933 -1.145 -0.859 -0.607 -0.261 -0.548 -0.519 -0.423 -0.260 -0.521 -0.161 -0.639 -0.264 0.223 -0.199] ;
    [-0.788 -0.658 -0.790 -0.669 -0.524 -0.673 -0.884 -0.670 -0.371 -0.230 -0.463 -0.423 -0.367 -0.253 -0.344 0.160 -0.455 -0.114 0.271 -0.018] ;
    [-0.835 -0.720 -0.807 -0.778 -0.729 -0.642 -0.997 -0.687 -0.323 0.033 -0.342 -0.260 -0.253 0.054 0.022 0.179 -0.290 -0.042 0.334 -0.035] ;
    [-0.616 -0.409 -0.482 -0.402 -0.291 -0.298 -0.613 -0.631 -0.235 -0.097 -0.382 -0.521 -0.344 0.022 0.179 0.634 -0.664 -0.584 -0.176 0.189] ;
    [-0.179 -0.209 -0.419 -0.439 -0.366 -0.335 -0.624 -0.453 -0.039 0.443 -0.192 -0.161 0.160 0.179 0.634 0.933 -0.324 -0.374 -0.057 0.257] ;
    [-1.499 -1.252 -1.330 -1.234 -1.176 -1.118 -1.383 -1.222 -0.646 -0.325 -0.720 -0.639 -0.455 -0.290 -0.664 -0.324 -1.078 -0.307 0.388 -0.346] ;
    [-0.771 -0.611 -0.805 -0.854 -0.758 -0.664 -0.912 -0.745 -0.327 -0.050 -0.247 -0.264 -0.114 -0.042 -0.584 -0.374 -0.307 0.200 0.815 -0.023] ;
    [-0.112 -0.146 -0.270 -0.253 -0.222 -0.200 -0.391 -0.349 0.196 0.589 0.155 0.223 0.271 0.334 -0.176 -0.057 0.388 0.815 1.339 0.661] ;
    [-1.196 -0.788 -1.076 -0.991 -0.771 -0.886 -1.278 -1.067 -0.374 -0.042 -0.222 -0.199 -0.018 -0.035 0.189 0.257 -0.346 -0.023 0.661 0.129]
]

const Full2aux = [
    [-5.44 -4.99 -5.80 -5.50 -5.83 -4.96 -4.95 -4.16 -3.57 -3.16 -3.11 -2.86 -2.59 -2.85 -2.41 -2.27 -3.60 -2.57 -1.95 -3.07] ; 
    [0 -5.46 -6.56 -6.02 -6.41 -5.32 -5.55 -4.91 -3.94 -3.39 -3.51 -3.03 -2.95 -3.30 -2.57 -2.89 -3.98 -3.12 -2.48 -3.45] ; 
    [0 0 -7.26 -6.84 -7.28 -6.29 -6.16 -5.66 -4.81 -4.13 -4.28 -4.02 -3.75 -4.10 -3.48 -3.56 -4.77 -3.98 -3.36 -4.25] ; 
    [0 0 0 -6.54 -7.04 -6.05 -5.78 -5.25 -4.58 -3.78 -4.03 -3.52 -3.24 -3.67 -3.17 -3.27 -4.14 -3.63 -3.01 -3.76] ; 
    [0 0 0 0 -7.37 -6.48 -6.14 -5.67 -4.91 -4.16 -4.34 -3.92 -3.74 -4.04 -3.40 -3.59 -4.54 -4.03 -3.37 -4.20] ; 
    [0 0 0 0 0 -5.52 -5.18 -4.62 -4.04 -3.38 -3.46 -3.05 -2.83 -3.07 -2.48 -2.67 -3.58 -3.07 -2.49 -3.32] ; 
    [0 0 0 0 0 0 -5.06 -4.66 -3.82 -3.42 -3.22 -2.99 -3.07 -3.11 -2.84 -2.99 -3.98 -3.41 -2.69 -3.73] ; 
    [0 0 0 0 0 0 0 -4.17 -3.36 -3.01 -3.01 -2.78 -2.76 -2.97 -2.76 -2.79 -3.52 -3.16 -2.60 -3.19] ; 
    [0 0 0 0 0 0 0 0 -2.72 -2.31 -2.32 -2.01 -1.84 -1.89 -1.70 -1.51 -2.41 -1.83 -1.31 -2.03] ; 
    [0 0 0 0 0 0 0 0 0 -2.24 -2.08 -1.82 -1.74 -1.66 -1.59 -1.22 -2.15 -1.72 -1.15 -1.87] ; 
    [0 0 0 0 0 0 0 0 0 0 -2.12 -1.96 -1.88 -1.90 -1.80 -1.74 -2.42 -1.90 -1.31 -1.90] ; 
    [0 0 0 0 0 0 0 0 0 0 0 -1.67 -1.58 -1.49 -1.63 -1.48 -2.11 -1.62 -1.05 -1.57] ; 
    [0 0 0 0 0 0 0 0 0 0 0 0 -1.68 -1.71 -1.68 -1.51 -2.08 -1.64 -1.21 -1.53] ; 
    [0 0 0 0 0 0 0 0 0 0 0 0 0 -1.54 -1.46 -1.42 -1.98 -1.80 -1.29 -1.73] ; 
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1.21 -1.02 -2.32 -2.29 -1.68 -1.33] ; 
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -0.91 -2.15 -2.27 -1.80 -1.26] ; 
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3.05 -2.16 -1.35 -2.25] ; 
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1.55 -0.59 -1.70] ; 
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -0.12 -0.97] ; 
    [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1.75] 
]

const Full2intMatrix = Symmetric(Full2aux)

# Declare an interaction model structure, includes the name, interaction matrix and a dictionary.
struct PF_model
    pf_name :: PFmodelname 
    int_matrix :: Array{Float64,2}
    dict :: Dict{Amin, Int8}
end


# Define the six available models.
const HP1_model = PF_model(HP1,HP1intMatrix,DictHP)
const HP2_model = PF_model(HP2,HP2intMatrix,DictHP)
const HP3_model = PF_model(HP3,HP3intMatrix,DictHP)
const HPNX_model = PF_model(HPNX,HPNXintMatrix,DictHP)
const hHPNX_model = PF_model(hHPNX,hHPNXintMatrix,DictHP)
const YhHX_model = PF_model(YhHX,YhHXintMatrix,DictYhHX)
const Full1_model = PF_model(Full1,Full1intMatrix,DictFull1)
const Full2_model = PF_model(Full2,Full2intMatrix,DictFull1)


# Deifine a protein structure.
struct Protein
    edo :: Matrix{Int16} # This encodes the amino acids' positions within the array.
    HPlist :: Array{Amin,1} # This encodes the protein's amino acid sequence.
    geometry :: geometries # This defines the type of geometry that the protein is embedded in.
end



@enum directions::Int8 begin
        forwards = 1
        backwards = 2
        nonetaken = 3    
end



const translate_HP = Dict{Amin,Amin}(CYS => P, 
MET => H,
PHE => H, 
ILE => H,
LEU => H,
VAL => H,
TRP => H,
TYR => P,
ALA => H,
GLY => H, 
THR => P,
SER => P, 
ASN => P, 
GLN => P,
ASP => P,
GLU => P,
HIS => P, 
ARG => P,
LYS => P,
PRO => H)
const translate_HPNX = Dict{Amin,Amin}(CYS => X, 
MET => H,
PHE => H, 
ILE => H,
LEU => H,
VAL => H,
TRP => H,
TYR => X,
ALA => H,
GLY => H, 
THR => X,
SER => X, 
ASN => X, 
GLN => X,
ASP => N,
GLU => N,
HIS => P, 
ARG => P,
LYS => P,
PRO => H)  
const translate_hHPNX = Dict{Amin,Amin}(CYS => X, 
MET => H,
PHE => H, 
ILE => H,
LEU => H,
VAL => h,
TRP => H,
TYR => X,
ALA => h,
GLY => H, 
THR => X,
SER => X, 
ASN => X, 
GLN => X,
ASP => N,
GLU => N,
HIS => P, 
ARG => P,
LYS => P,
PRO => H)
const translate_YhHX = Dict{Amin,Amin}(CYS => H, 
MET => H,
PHE => H, 
ILE => H,
LEU => H,
VAL => h,
TRP => X,
TYR => Y,
ALA => h,
GLY => Y, 
THR => X,
SER => Y, 
ASN => Y, 
GLN => X,
ASP => X,
GLU => Y,
HIS => Y, 
ARG => Y,
LYS => X,
PRO => X)


"""
    translate_HPlist(HPlist,pfmodel,translate)
Given an amino acid sequence `HPlist`, a model `pfmodel` and a boolean value `translate`; returns the equivalent list for the given model.
"""
function translate_HPlist(HPlist::Vector{Amin},pfmodel::PFmodelname,translate::Bool)
    if translate
        len = length(HPlist)
        HPlist_aux = fill(CYS,len)
        if (pfmodel == HP1) || (pfmodel == HP2) || (pfmodel == HP3)
            for k in 1:len
                HPlist_aux[k] = translate_HP[HPlist[k]]
            end
    
        elseif pfmodel == HPNX
            for k in 1:len
                HPlist_aux[k] = translate_HPNX[HPlist[k]]
            end
    
        elseif pfmodel == hHPNX
            for k in 1:len
                HPlist_aux[k] = translate_hHPNX[HPlist[k]]
            end
    
        elseif pfmodel == YhHX
            for k in 1:len
                HPlist_aux[k] = translate_YhHX[HPlist[k]]
            end
        end
        return HPlist_aux

    else
        return HPlist
    end
end











"""
    dirsf(vec)
Given a vector with integer entries; returns the equivalent `dirs` vector.
"""
function dirsf(vec::Array{Int8,1})
    l = length(vec)
    dirsVec = fill(nonetaken,l) # Declare the equivalent `dirs` vector.
    for k in 1:l
        val = vec[k]
        if val == 1
            dirsVec[k] = forwards
        elseif val == 2
            dirsVec[k] = backwards
        end
    end
    return dirsVec
end

















"""
    periodicInd(A,indices,dim)

Given a 3D array `A`, a couple/triad of indices representing coordinates `indices` and a dimension `dim`; returns the indices (coordinate) for 
the equivalent array with periodic boundary conditions.
"""
function periodicInd(A,indices,dim::Int) 
    if dim == 2
        lx,ly = size(A)
        ix,iy = indices
        Ix = mod1(ix,lx)
        Iy = mod1(iy,ly)
        return Int16[Ix,Iy]
    else
        lx,ly,lz = size(A)
        ix,iy,iz = indices
        Ix = mod1(ix,lx)
        Iy = mod1(iy,ly)
        Iz = mod1(iz,lz)
        return Int16[Ix,Iy,Iz]
    end
end










"""
    periodicArr(A,indices,dim)

Given a 3D array `A`, a couple/triad of indices representing coordinates `indices` and a dimension `dim`; returns the value of the position `indices` in 
the equivalent array with periodic boundary conditions.
"""
function periodicArr(A,indices,dim::Int)
    if dim == 2
        lx,ly = size(A)
        ix,iy = indices
        Ix = mod1(ix,lx)
        Iy = mod1(iy,ly)
        return A[Ix,Iy]
    else
        lx,ly,lz = size(A)
        ix,iy,iz = indices
        Ix = mod1(ix,lx)
        Iy = mod1(iy,ly)
        Iz = mod1(iz,lz)
        return A[Ix,Iy,Iz]
    end
end










"""
    makeLattice(N,edo,HPlist)

Given a value for the lattice size `N`, a matrix encoding the aminoacids positions `edo`, and an array containing 
the sequence of aminoacids `HPlist`; creates a 2D/3D array containing the amino acid sequence.
"""
function makeLattice(N::Int,edo,HPlist::Array{Amin,1})
    if length(edo[1,:]) == 3
        red = zeros(Int8,(N,N,N))
        for k in 1:length(HPlist)
            x,y,z = periodicInd(red,edo[k,:],3)
            red[x,y,z] = Int(HPlist[k])
        end    
        return red
    else
        red = zeros(Int8,(N,N))
        for k in 1:length(HPlist)
            x,y = periodicInd(red,edo[k,:],2)
            red[x,y] = Int(HPlist[k])
        end    
        return red
    end
end














# Next, generalize the concept of nearest neighbors, according to the geometry.
"""
    nearestNeighbors(red,inds,geometry)

Given a 2D/3D array `red`, a poition `inds`, and a geometry `geometry`; returns the values of the nearest 
neighbors to the given position. 
"""
function nearestNeighbors(red,inds,geometry::geometries)
    A = red

    if geometry == cubic
        x,y,z = inds
        nn = Int8[periodicArr(A,[x-1,y,z],3),periodicArr(A,[x,y+1,z],3),periodicArr(A,[x+1,y,z],3),periodicArr(A,[x,y-1,z],3),
        periodicArr(A,[x,y,z+1],3),periodicArr(A,[x,y,z-1],3)] # Last two neighbors fall outside x-y plane.
        return nn
    

    elseif geometry == fcc # Fcc geometry has 12 topological nearest neighbors.
        i,j,k = inds
        nn = Int8[periodicArr(A,[i+1,j,k],3),periodicArr(A,[i-1,j,k],3)
        ,periodicArr(A,[i,j+1,k],3),periodicArr(A,[i,j-1,k],3)
        ,periodicArr(A,[i,j,k+1],3),periodicArr(A,[i,j,k-1],3)
        ,periodicArr(A,[i+1,j-1,k],3),periodicArr(A,[i-1,j+1,k],3)
        ,periodicArr(A,[i,j+1,k-1],3),periodicArr(A,[i,j-1,k+1],3)
        ,periodicArr(A,[i-1,j,k+1],3),periodicArr(A,[i+1,j,k-1],3)]
        return nn
    

    elseif geometry == square2D
        x,y = inds
        nn = Int8[periodicArr(A,[x-1,y],2),periodicArr(A,[x,y+1],2),periodicArr(A,[x+1,y],2),periodicArr(A,[x,y-1],2)]
        return nn


    elseif geometry == triangular2D
        x,y = inds
        nn = Int8[periodicArr(A,[x-1,y],2),periodicArr(A,[x,y+1],2)
        ,periodicArr(A,[x+1,y+1],2),periodicArr(A,[x+1,y],2)
        ,periodicArr(A,[x,y-1],2),periodicArr(A,[x-1,y-1],2)]
        return nn
    end

end















"""
    nearestNeighborsCoords(red,inds,geometry)

Given a 2D/3D array `red`,a position `inds`, and a geometry `geometry`; returns the coordinates of the nearest neighbors 
to the given position. 
"""
function nearestNeighborsCoords(red,inds,geometry::geometries)
    A = red

    if geometry == cubic
        x,y,z = inds 
        nnc = Vector{Int16}[periodicInd(A,[x-1,y,z],3),periodicInd(A,[x,y+1,z],3),periodicInd(A,[x+1,y,z],3),periodicInd(A,[x,y-1,z],3),
        periodicInd(A,[x,y,z+1],3),periodicInd(A,[x,y,z-1],3)] # Last two neighbors fall outside x-y plane.
        return nnc
    
    
    elseif geometry == fcc # First six neighbors are in the same x-y plane. The remaining six are outside.
        i,j,k = inds
        nnc = Vector{Int16}[periodicInd(A,[i+1,j,k],3),periodicInd(A,[i-1,j,k],3)
        ,periodicInd(A,[i,j+1,k],3),periodicInd(A,[i,j-1,k],3)
        ,periodicInd(A,[i,j,k+1],3),periodicInd(A,[i,j,k-1],3)
        ,periodicInd(A,[i+1,j-1,k],3),periodicInd(A,[i-1,j+1,k],3)
        ,periodicInd(A,[i,j+1,k-1],3),periodicInd(A,[i,j-1,k+1],3)
        ,periodicInd(A,[i-1,j,k+1],3),periodicInd(A,[i+1,j,k-1],3)]
        return nnc
    

    elseif geometry == square2D
        x,y = inds
        nnc  = Vector{Int16}[periodicInd(A,[x-1,y],2),periodicInd(A,[x,y+1],2),periodicInd(A,[x+1,y],2),periodicInd(A,[x,y-1],2)]
        return nnc


    elseif geometry == triangular2D
        x,y = inds
        nnc = Vector{Int16}[periodicInd(A,[x-1,y],2),periodicInd(A,[x,y+1],2)
        ,periodicInd(A,[x+1,y+1],2),periodicInd(A,[x+1,y],2)
        ,periodicInd(A,[x,y-1],2),periodicInd(A,[x-1,y-1],2)]
        return nnc
    end
end

















"""
    sharedNeighborsCoords(red,inds,indsp,geometry)

Given a 2D/3D array `red`, a couple of coordinates `inds` and `indsp`, and a geometry; returns the coordinates for the empty shared topological
neighbors of `inds,indsp` for the given geometry. 
"""
function sharedNeighborsCoords(red,inds,indsp,geometry::geometries)
    nnc = nearestNeighborsCoords(red,inds,geometry) # Coordiantes of nearest neigbors to our coord `ind`.
    nn = nearestNeighbors(red,inds,geometry) # Value of nearest neighbors to our index.
    nncp = nearestNeighborsCoords(red,indsp,geometry) # Coordiantes of nearest neigbors to our coord `indsp`.
    

    sharedNC = Vector{Int16}[] # Empty shared neighbor spaces will have a value of zero.
    for i in 1:length(nnc)
        if (nnc[i] ∈ nncp) && (nn[i] == 0)
            push!(sharedNC,nnc[i]) # These are the coordinates where we might move the index corresponding to `inds`
        end
    end

    return sharedNC
end
















"""
    excludedNeighborsCoords(red,inds,indsp,geometry)

Given a 2D/3D array `red`, a couple of coordinates `inds` and `indsp`, and a geometry; returns the coordinates for the empty not-shared neighbor spaces by 
`inds,indsp` for the given geometry.  
"""
function excludedNeighborsCoords(red,inds,indsp,geometry::geometries)
    nnc = nearestNeighborsCoords(red,inds,geometry) # Coordiantes of nearest neigbors to our index `ind`.
    nn = nearestNeighbors(red,inds,geometry) # Value of nearest neighbors to our index.
    nncp = nearestNeighborsCoords(red,indsp,geometry) # Coordiantes of nearest neigbors to our index `indsp`.
    

    notsharedNC = Vector{Int16}[] # Empty not shared neighbor spaces will have a value of zero.
    for i in 1:length(nnc)
        if (nnc[i] ∉ nncp) && (nn[i] == 0)
            push!(notsharedNC,nnc[i]) # These are the coordinates where we might move the index corresponding to `inds`
        end
    end

    return notsharedNC
end























# Now I write a function which checks for the validity of a configuration, which also depends on the geometry of the configuration.
"""
    validConf(N,ind,edo,HPlist,dir,geometry)

Given a 2D/3D array of size `N`, a matrix encoding the aminoacids positions `edo`, an index on the chain `ind` , an array containing 
the sequence of aminoacids `HPlist`, a direction `dir`, and a geometry; determines whether the protein structure is valid or not.
"""
function validConf(N::Int,ind,edo,HPlist::Vector{Amin},dir::directions,geometry::geometries)
    # I set up the protein within the array.
    red = makeLattice(N,edo,HPlist)

    # Next, I iterate over the protein´s vertices, checking wether each pair of positions is in each other´s list of nearest neighbors.
    # The direction in which I check the structure is determined by `dir`.
    
    ans = true

    
    if (geometry == fcc) || (geometry == cubic) 
        if dir == backwards # Check backwards.
            if ind != 1
                for j in ind:-1:2
                    x1,y1,z1 = edo[j,:]
                    x2,y2,z2 = periodicInd(red,edo[j-1,:],3) 
                    nnc = nearestNeighborsCoords(red,[x1,y1,z1],geometry) # Nearest neigbors to the position `[x1,y1,z1]`.
                    
                    if ([x2,y2,z2] ∉ nnc)  # Configuration is not valid.
                        ans = false
                        break
                    end
                end
            end
    
    
        elseif dir == forwards # Check forwards.
            if ind != length(HPlist)
                for j in ind:length(HPlist)-1
                    x1,y1,z1 = edo[j,:]
                    x2,y2,z2 = periodicInd(red,edo[j+1,:],3) 
                    nnc = nearestNeighborsCoords(red,[x1,y1,z1],geometry) # Nearest neigbors to the position `[x1,y1,z1]`.
                    

                    if ([x2,y2,z2] ∉ nnc)  # Configuration is not valid.
                        ans = false
                        break
                    end
                end
            end

        end

    elseif (geometry == square2D) || (geometry == triangular2D) 
        if dir == backwards
            if ind != 1
                for j in ind:-1:2
                    x1,y1 = edo[j,:] 
                    x2,y2 = periodicInd(red,edo[j-1,:],2) 
                    nnc = nearestNeighborsCoords(red,[x1,y1],geometry) # Nearest neigbors to the position `[x1,y1]`.
                    
                    if ([x2,y2] ∉ nnc) # Configuration is not valid.
                        ans = false
                        break
                    end
                end
            end
    
    
        elseif dir == forwards
            if ind != length(HPlist)
                for j in ind:length(HPlist)-1
                    x1,y1 = edo[j,:]
                    x2,y2 = periodicInd(red,edo[j+1,:],2) # Make the indices periodic.
                    nnc = nearestNeighborsCoords(red,[x1,y1],geometry) # Nearest neigbors to the position `[x1,y1]`.
                    

                    if ([x2,y2] ∉ nnc)  # Configuration is not valid.
                        ans = false
                        break
                    end
                end
            end

        end

    end


    return ans
end


















"""
    countFirst(N,edo,HPlist,geometry)

Given a 2D/3D array size `N`, a matrix encoding the aminoacids positions `edo`,an array containing the sequence of aminoacids `HPlist`, and
a geometry; counts the number of possible pull moves for the first amino acid, recording the possible positions to which we might move the first
(and second in the case of the square geometry) aminoacid to.
"""
function countFirst(N::Int,edo,HPlist::Vector{Amin},geometry::geometries)  
    
    if geometry == square2D || geometry == cubic
        red = makeLattice(N,edo,HPlist)
        ind = 1
        # Find out which spaces are free for `ind+1` to move into.
        coordinatesp = excludedNeighborsCoords(red,edo[ind,:],edo[ind+1,:],geometry)

        # Declare the arrays that will contain the coordinates for each of the different pull-moves.
        coordinates1 = Vector{Int16}[] # These are the coordinates for `ind+1`
        coordinates2 = Vector{Int16}[] # These are the coordinates for `ind`
        
        for coord1 in coordinatesp # This are the coordinates to which I should be able to move monomer `ind+1` into.
            coordsaux = excludedNeighborsCoords(red,coord1,edo[ind+1,:],geometry)
            l = length(coordsaux)
            append!(coordinates1,Vector{Int16}[coord1 for i in 1:l])
            append!(coordinates2,coordsaux)
        end
    
        numpull = length(coordinates1)
        # Now we have the possible coordinates for the first two monomers, as well as the number of possible pull moves for the first 
        # amino acid. Fo easier use, I turn the arrays `coordinates1,coordinates2` into matrices.
    
        coordinates1p = transpose(hcat(coordinates1...))
        coordinates2p = transpose(hcat(coordinates2...))
    
        return (numpull,coordinates1p,coordinates2p)
        
    
    elseif (geometry == fcc) || (geometry == triangular2D)
        red = makeLattice(N,edo,HPlist)
        ind = 1

        # Find out which spaces are free for `ind` to move into.
        coordinates1 = excludedNeighborsCoords(red,edo[ind,:],edo[ind+1,:],geometry)

        numpull = length(coordinates1) # Compute the number of pull-moves.
        # Now we have the possible coordinates for the first monomer, as well as the number of possible pull moves for the first 
        # amino acid. Fo easier use, I turn the array into a matrix.

        coordinates1p = transpose(hcat(coordinates1...))
        return (numpull,coordinates1p)
    end
end


















"""
    countLast(N,edo,HPlist,geometry) 

Given a 2D/3D array size `N`, a matrix encoding the aminoacids positions `edo`,an array containing the sequence of aminoacids `HPlist`, and
a geometry; counts the number of possible pull moves for the last amino acid in the chain, recording the possible positions.
"""
function countLast(N::Int,edo,HPlist::Vector{Amin},geometry::geometries) 

    if geometry == square2D || geometry == cubic
        red = makeLattice(N,edo,HPlist)
        ind = length(HPlist)
        # Find out which spaces are free for `ind-1` to move into.
        coordinatesp = excludedNeighborsCoords(red,edo[ind,:],edo[ind-1,:],geometry)

        # Declare the arrays that will contain the coordinates for each of the different pull-moves.
        coordinates1 = Vector{Int16}[] # These are the coordinates for `ind-1`
        coordinates2 = Vector{Int16}[] # These are the coordinates for `ind`
        
        for coord1 in coordinatesp # This are the coordinates to which I should be able to move monomer `ind-1` into.
            coordsaux = excludedNeighborsCoords(red,coord1,edo[ind-1,:],geometry)
            l = length(coordsaux)
            append!(coordinates1,Vector{Int16}[coord1 for i in 1:l])
            append!(coordinates2,coordsaux)
        end
    
        numpull = length(coordinates1)
        # Now we have the possible coordinates for the last two monomers, as well as the number of possible pull moves for the last 
        # amino acid. Fo easier use, I turn the arrays `coordinates1,coordinates2` into matrices.
    
        coordinates1p = transpose(hcat(coordinates1...))
        coordinates2p = transpose(hcat(coordinates2...))
        
        return (numpull,coordinates1p,coordinates2p)

    
    elseif (geometry == fcc) || (geometry == triangular2D)
        red = makeLattice(N,edo,HPlist)
        ind = length(HPlist)

        # Find out which spaces are free for `ind` to move into.
        coordinates1 = excludedNeighborsCoords(red,edo[ind,:],edo[ind-1,:],geometry)

        numpull = length(coordinates1) # Compute the number of pull-moves.
        # Now we have the possible coordinates for the last monomer, as well as the number of possible pull moves for the last
        # amino acid. Fo easier use, I turn the array into a matrix.

        coordinates1p = transpose(hcat(coordinates1...))
        return (numpull,coordinates1p)
    end
end

















"""
    countMiddle(N,ind,edo,HPlist,geometry) 

Given a 2D/3D array size `N`, an index on the chain `ind`, a matrix encoding the aminoacids positions `edo`, an array containing 
the sequence of aminoacids `HPlist`, and a geometry; counts the number of possible pull moves for the `ìnd`-th aminoacid, and 
stores the possible coordinates.
"""
function countMiddle(N::Int,ind,edo,HPlist::Vector{Amin},geometry::geometries)
    
    if (geometry == fcc) || (geometry == triangular2D)
        # I set up the protein within the array.
        red = makeLattice(N,edo,HPlist)
    
        # Now we check whether a pull-move is possible for the selected vertex. A shared neighbor space between `ind` and either `ind+1` or `ind-1`
        # must be empty.

        # Next we store the possible free spaces for `ind` when pulling backwards (as in moving only the previous monomers).
        coordinates1 = sharedNeighborsCoords(red,edo[ind,:],edo[ind+1,:],geometry)
    
        # Store possible free spaces when pulling forwards.
        coordinatesi = sharedNeighborsCoords(red,edo[ind,:],edo[ind-1,:],geometry)
    
    
        # Now I have in `coordinates1,coordinatesi` the potential new positions for `ind`. The length of these arrays is 
        # the number of possible pull moves por the given index.
        numpull = length(coordinates1)+length(coordinatesi)

        # I turn the 1 dimensional arrays into matrices for easier use.
        coordinates1p = transpose(hcat(coordinates1...))
        coordinatesip = transpose(hcat(coordinatesi...))

        return (numpull,coordinates1p,coordinatesip) # Returns everything neccesary to reproduce the final configuration.


    elseif geometry == square2D
        # Set up the protein within the array.
        red = makeLattice(N,edo,HPlist)

        # Next, we store the possible coordinates to which we might be able to move `ind-1` or `ind+1`.
        coordinatesaux1 = Vector{Int16}[]
        coordinatesaux2 = Vector{Int16}[]
        for coord in nearestNeighborsCoords(red,edo[ind,:],geometry) # Leaving an extra space here. !!!!!!!!!
            if (coord == edo[ind-1,:]) || (red[CartesianIndex(coord[1],coord[2])] == 0)
                push!(coordinatesaux1,coord)
            end

            if (coord == edo[ind+1,:]) || (red[CartesianIndex(coord[1],coord[2])] == 0)
                push!(coordinatesaux2,coord)
            end
        end

        # First consider the case when we are pulling the protein backwards (as in pulling only the preceeding monomers to our index `ind`).
        coordinates1 = Vector{Int16}[] # Coordinates for `ind`
        coordinates2 = Vector{Int16}[] # Coordinates for `ind-1`.

        for coord2 in coordinatesaux1
            coord1 = sharedNeighborsCoords(red,coord2,edo[ind+1,:],geometry)
            append!(coordinates1,coord1)
            append!(coordinates2,Vector{Int16}[coord2 for i in 1:length(coord1)])
        end

        # Now consider the case when  we are pulling the chain forwards.
        coordinatesi = Vector{Int16}[] # Coordinates for `ind`
        coordinatesii = Vector{Int16}[] # Coordinates for `ind+1`.

        for coordii in coordinatesaux2
            coordi = sharedNeighborsCoords(red,coordii,edo[ind-1,:],geometry)
            append!(coordinatesi,coordi)
            append!(coordinatesii,Vector{Int16}[coordii for i in 1:length(coordi)])
        end

        numpull = length(coordinates1)+length(coordinatesi)

        # I turn the 1 dimensional arrays into matrices for easier use.
        coordinates1p = transpose(hcat(coordinates1...))
        coordinates2p = transpose(hcat(coordinates2...))
        coordinatesip = transpose(hcat(coordinatesi...))
        coordinatesiip = transpose(hcat(coordinatesii...))

        return (numpull,coordinates1p,coordinates2p,coordinatesip,coordinatesiip) # Returns everything neccesary to reproduce the final configuration.
    end

end




















"""
    countpull(N,edo,HPlist,geometry)

Given a 2D/3D array size `N`, a matrix encoding the aminoacids positions `edo`, an array containing 
the sequence of aminoacids `HPlist`, and a geometry; counts all of the possible pull moves. It also outputs the necessary
coordinates to perform the listed moves.
"""
function countpull(N::Int,edo,HPlist::Vector{Amin},geometry::geometries)

    if (geometry == fcc) || (geometry == triangular2D)
        # Count the moves and store the coordinates for the first monomer.
        t1,coords1 = countFirst(N,edo,HPlist,geometry)

        # Count the moves and store the coordinates for the last monomer.
        t2,coords3 = countLast(N,edo,HPlist,geometry)

        # Count the moves and store the coordinates for the middle amino acids. 
        t3 = 0
        tmid = 0
        indexb = Int16[] # Stores the value of the middle indices that have non empty arrays of possible backwards pull moves.
        indexf = Int16[] # Stores the value of the middle indices that have non empty arrays of possible forwards pull moves.
        npullpindexb = Int16[] # Stores the number of backwards pull moves for each middle index.
        # Length may be shorter than the number of middle indices given that I only store the coordinates for non empty arrays.
        npullpindexf = Int16[] # Stores the number of forwards pull moves for each middle index.
        middlecoords1 = Matrix{Int16}[] # Contains the coordinates for `ind`, where `ind` is an index from the middle of the chain.
        middlecoords3 = Matrix{Int16}[] # Contains the coordinates for `ind` when the chain is pulled forwards.

        for j in 2:length(HPlist)-1
            tj,coords5,coords7 = countMiddle(N,j,edo,HPlist,geometry)
            # I have all the neccesary information. But first I need to check if the given index
            # posseses at least a possible pull move.
            t3 = t3+tj
            if isempty(coords5) == false
                tmid = tmid+size(coords5)[1]
                push!(indexb,j)
                push!(npullpindexb,size(coords5)[1])
                push!(middlecoords1,coords5)
            end
    
            if isempty(coords7) == false
                push!(indexf,j)
                push!(npullpindexf,size(coords7)[1])
                push!(middlecoords3,coords7)
            end
        end

        totalpull = t1+t2+t3
        matrixb = hcat(indexb,npullpindexb)
        matrixf = hcat(indexf,npullpindexf)

        return (totalpull,tmid,matrixb,matrixf,coords1,coords3,middlecoords1,middlecoords3)

    elseif  geometry == square2D
        # Count the moves and store the coordinates for the first monomer.
        t1,coords2,coords1 = countFirst(N,edo,HPlist,geometry)

        # Count the moves and store the coordinates for the last monomer.
        t2,coords4,coords3 = countLast(N,edo,HPlist,geometry)

        # Count the moves and store the coordinates for the middle amino acids. 
        t3 = 0
        tmid = 0
        indexb = Int16[] # Stores the value of the middle indices that have non empty arrays of possible backwards pull moves.
        indexf = Int16[] # Stores the value of the middle indices that have non empty arrays of possible forwards pull moves.
        npullpindexb = Int16[] # Stores the number of backwards pull moves for each middle index.
        # Length may be shorter than the number of middle indices given that I only store the coordinates for non empty arrays.
        npullpindexf = Int16[] # Stores the number of forwards pull moves for each middle index.
        middlecoords1 = Matrix{Int16}[] # Contains the coordinates for `ind`, where `ind` is an index from the middle of the chain.
        middlecoords2 = Matrix{Int16}[] # Contains the coordinates for `ind-1`.
        middlecoords3 = Matrix{Int16}[] # Contains the coordinates for `ind`.
        middlecoords4 = Matrix{Int16}[] # Contains the coordinates for `ind+1`.
        for j in 2:length(HPlist)-1
            tj,coords5,coords6,coords7,coords8 = countMiddle(N,j,edo,HPlist,geometry)
            # I have all the neccesary information. But first I need to check if the given index
            # posseses at least a possible pull move.
            t3 = t3+tj
            if isempty(coords5) == false
                tmid = tmid+size(coords5)[1]
                push!(indexb,j)
                push!(npullpindexb,size(coords5)[1])
                push!(middlecoords1,coords5)
                push!(middlecoords2,coords6)
            end
    
            if isempty(coords7) == false
                push!(indexf,j)
                push!(npullpindexf,size(coords7)[1])
                push!(middlecoords3,coords7)
                push!(middlecoords4,coords8)
            end
        end

        totalpull = t1+t2+t3

        matrixb = hcat(indexb,npullpindexb)
        matrixf = hcat(indexf,npullpindexf)

        return (totalpull,tmid,matrixb,matrixf,coords1,coords2,coords3,coords4,middlecoords1,middlecoords2,middlecoords3,middlecoords4)
    end

end



















# Before writing the function that performs the actual move, I need a function which returns the index being pullled
# in one of the `middlecoords`matrices from my function `countpull3D`.
"""
    middleInd(mp,matrix)

Given a number `mp`, and a matrix `matrix` containing the values of the middle indices with possible pull moves in 
the first column and the number of pull moves for said indices in the second column;returns the corresponding 
index `ind` being moved and the position in one of my `middlecoords[ind][position,:]` matrices.
"""
function middleInd(mp,matrix)
    indexd = matrix[:,1] # This array contains the value of the indices with viable pull moves.
    npullpindex = matrix[:,2] # This array contains the number of pull moves for each index in `indexd`.
    l = length(npullpindex)
    nb = ones(Int16,l) # `nb` will store the number of accumulated pull moves for each of the indices.
    nb[1] = npullpindex[1]
    for i in 2:l # Fill `nb`.
        nb[i] = npullpindex[i]+nb[i-1]
    end

    indm = 0 # This is the index being pulled.
    ind = 0 # This is the index for the `middlecoords[ind][pos,:]` matrices.
    pos = 0 # This is the position of the right coordinates.
    
    
    for i in 1:l
        if mp ≤ nb[i]
            indm = indexd[i]
            ind = i
            break
        end
    end

    if mp ≤ nb[1]
        pos = mp
    else
        pos = mp-nb[ind-1]
    end
    
    return (indm,ind,pos)
end




















"""
    pullMove(N,edo,HPlist,geometry)  

Given a 2D/3D array size `N`, a matrix encoding the aminoacids positions `edo`, an array containing 
the sequence of aminoacids `HPlist`, and a geometry; chooses an index and performs a full pull move for the generated index. 
The function also returns the final state `newedo`, a number of possible pull moves `totalpull` for the 
initial configuration, the position of the pulled monomer `indpulled`, and the new position for the pulled index `newcoord`.
"""
function pullMove(N::Int,edo,HPlist::Vector{Amin},geometry::geometries)   

    if geometry == fcc || geometry == triangular2D
        # `newedo` will contain the amino acids´ final positions.
        newedo = copy(edo)

        # Generate the list of possible pull moves.
        totalpull,tmid,matrixb,matrixf,coords1,coords3,middlecoords1,middlecoords3 = countpull(N,edo,HPlist,geometry)

        # Randomly choose one of the pull moves.
        m = rand(1:totalpull)
    
        # Make subdivisions according to the type of move.
        if isempty(coords1)
            s1 = 0
        else
            s1 = size(coords1)[1] 
        end

        if isempty(coords3)
            s2 = 0
            s3 = s1
    
        else
            s2 = s1+1
            s3 = s1+(size(coords3)[1])
        end

        s4 = s3+1
        s5 = s3+tmid
        s6 = s5+1

        # Type of pull move depends on the value of `m`.
        if m ≤ s1 # First monomer is moved.
            ind = 1
            newedo[ind,:] = coords1[m,:] # Record the change.
            indpulled = ind # Function returns the index pulled in the move.
            newcoord = newedo[ind,:] # Function returns the new coordinate for the pulled index.
            dirpulled = forwards
    
            # If all went well, a move has been done, now I have to check whether the current configuration is valid.
            stateconf = validConf(N,ind,newedo,HPlist,forwards,geometry)
            if stateconf == false
                for k in (ind+1):length(HPlist)
                    if stateconf == false
                        newedo[k,:] = edo[k-1,:] # New position. 
                        stateconf = validConf(N,k,newedo,HPlist,forwards,geometry) # Check whether the new configuration is valid.
                    else
                        break
                    end
                end
            end
    
    
    
    
        elseif  s2 ≤ m ≤ s3
            mp = m-(s2-1) # `mp is the position` for the moves corresponding to the final monomer in the chain.
            ind = length(HPlist)
            newedo[ind,:] = coords3[mp,:] # Record the change.
            indpulled = ind
            newcoord = newedo[ind,:]
            dirpulled = backwards
    
            # If all went well, a move has been done, now I have to check whether the current configuration is valid.
            stateconf = validConf(N,ind,newedo,HPlist,backwards,geometry)
            if stateconf == false
                for k in (ind-1):-1:1
                    if stateconf == false
                        newedo[k,:] = edo[k+1,:] # New position. 
                        stateconf = validConf(N,k,newedo,HPlist,backwards,geometry) # Check whether the new configuration is valid.
                    else
                        break
                    end
                end
            end
            
    
    
    
        elseif s4 ≤ m ≤ s5
            mp = m-(s4-1)
            indm,ind,pos = middleInd(mp,matrixb)
            newedo[indm,:] = middlecoords1[ind][pos,:]
            indpulled = indm
            newcoord = newedo[indm,:]
            dirpulled = backwards
    
            # If all went well, a move has been done, now I have to check whether the current configuration is valid.
            stateconf = validConf(N,indm,newedo,HPlist,backwards,geometry)
            if stateconf == false
                for k in (indm-1):-1:1
                    if stateconf == false
                        newedo[k,:] = edo[k+1,:] # New position.
                        stateconf = validConf(N,k,newedo,HPlist,backwards,geometry) # Check whether the new configuration is valid.
                    else
                        break
                    end
                end
            end
            
    
    
        elseif m ≥ s6 
            mp = m-(s6-1)
            indm,ind,pos = middleInd(mp,matrixf)
            newedo[indm,:] = middlecoords3[ind][pos,:]
            indpulled = indm
            newcoord = newedo[indm,:]
            dirpulled = forwards
    
            # If all went well, a move has been performed, now I have to check whether the current configuration is valid.
            stateconf = validConf(N,indm,newedo,HPlist,forwards,geometry)
            if stateconf ==  false
                for k in (indm+1):length(HPlist)
                    if stateconf == false
                        newedo[k,:] = edo[k-1,:] # New position. 
                        stateconf = validConf(N,k,newedo,HPlist,forwards,geometry) # Check whether the new configuration is valid.
                    else
                        break
                    end
                end
            end
    
    
    
        end
    
        return (newedo,totalpull,indpulled,newcoord,dirpulled) # Returns everything neccesary to reproduce the final configuration and implement
        # the Metropolis algorithm.


    elseif geometry == square2D
        # `newedo` will contain the amino acids´ final positions.
        newedo = copy(edo)

        # Generate the list of possible pull moves.
        totalpull,tmid,matrixb,matrixf,coords1,coords2,coords3,coords4,middlecoords1,middlecoords2,middlecoords3,middlecoords4 = countpull(N,edo,HPlist,geometry)

        # Randomly choose one of the pull moves.
        m = rand(1:totalpull)
    
        # Make subdivisions according to the type of move.
        if isempty(coords1)
            s1 = 0
        else
            s1 = size(coords1)[1] 
        end

        if isempty(coords3)
            s2 = 0
            s3 = s1
    
        else
            s2 = s1+1
            s3 = s1+(size(coords3)[1])
        end

        s4 = s3+1
        s5 = s3+tmid
        s6 = s5+1
    

        # Type of pull move depends on the value of `m`.
        if m ≤ s1 # First monomer is moved.
            ind = 1
            newedo[ind,:] = coords1[m,:] # Record the change.
            newedo[ind+1,:] = coords2[m,:] # Record the change.
            indpulled = ind # Function returns the index pulled in the move.
            newcoord1 = newedo[ind,:] # Function returns the new coordinate for the pulled index.
            newcoord2 = newedo[ind+1,:] 
            dirpulled = forwards
    
    
            # Now I have to check whether the current configuration is valid.
            stateconf = validConf(N,ind+1,newedo,HPlist,forwards,geometry)
            for k in ind+2:length(HPlist)
                if stateconf == false
                    newedo[k,:] = edo[k-2,:] # New position. 
                    stateconf = validConf(N,k,newedo,HPlist,forwards,geometry) # Check whether the new configuration is valid.
                else
                    break
                end
            end
            
    
    
        elseif  s2 ≤ m ≤ s3
            mp = m-(s2-1) # `mp is the position` for the moves corresponding to the final monomer in the chain.
            ind = length(HPlist)
            newedo[ind,:] = coords3[mp,:] # Record the change.
            newedo[ind-1,:] = coords4[mp,:] # Record the change.
            indpulled = ind # Function returns the index pulled in the move.
            newcoord1 = newedo[ind,:] # Function returns the new coordinate for the pulled index.
            newcoord2 = newedo[ind-1,:] 
            dirpulled = backwards
    
            # Now I have to check whether the current configuration is valid.
            stateconf = validConf(N,ind-1,newedo,HPlist,backwards,geometry) 
            for k in ind:-1:3
                if stateconf == false
                    newedo[k-2,:] = edo[k,:] 
                    stateconf = validConf(N,k-2,newedo,HPlist,backwards,geometry) # Check whether the new configuration is valid.
                else
                    break
                end
            end
    
    
        elseif s4 ≤ m ≤ s5
            mp = m-(s4-1)
            indm,ind,pos = middleInd(mp,matrixb)
            newedo[indm,:] = middlecoords1[ind][pos,:]
            indpulled = indm # Function returns the index pulled in the move.
            newcoord1 = newedo[indm,:] # Function returns the new coordinate for the pulled index. 
            dirpulled = backwards
    
            if validConf(N,indm,newedo,HPlist,backwards,geometry) != true
                newedo[indm-1,:] = middlecoords2[ind][pos,:] # Record the change.
            end
            newcoord2 = newedo[indm-1,:]
    
            # If all went well, a move has been done, now I have to check whether the current configuration is valid.
            stateconf = validConf(N,indm-1,newedo,HPlist,backwards,geometry)
            if stateconf == false
                for k in (indm-2):-1:1
                    if stateconf == false
                        newedo[k,:] = edo[k+2,:] # New position.
                        stateconf = validConf(N,k,newedo,HPlist,backwards,geometry) # Check whether the new configuration is valid.
                    else
                        break
                    end
                end
            end
            
    
    
        elseif m ≥ s6 
            mp = m-(s6-1)
            indm,ind,pos = middleInd(mp,matrixf)
            newedo[indm,:] = middlecoords3[ind][pos,:]
            indpulled = indm # Function returns the index pulled in the move.
            newcoord1 = newedo[indm,:] # Function returns the new coordinate for the pulled index. 
            dirpulled = forwards
            
    
            if validConf(N,indm,newedo,HPlist,forwards,geometry) != true
                newedo[indm+1,:] = middlecoords4[ind][pos,:] # Record the change.
            end
            newcoord2 = newedo[indm+1,:]
    
    
            # If all went well, a move has been performed, now I have to check whether the current configuration is valid.
            stateconf = validConf(N,indm+1,newedo,HPlist,forwards,geometry)
            if stateconf ==  false
                for k in (indm+2):length(HPlist)
                    if stateconf == false
                        newedo[k,:] = edo[k-2,:] # New position. 
                        stateconf = validConf(N,k,newedo,HPlist,forwards,geometry) # Check whether the new configuration is valid.
                    else
                        break
                    end
                end
            end
    
    
    
        end
    

        return (newedo,totalpull,indpulled,newcoord1,newcoord2,dirpulled) # Returns everything neccesary to reproduce the final configuration and implement
        # the Metropolis algorithm.

    end
        
    
end


















"""
    reconstructStates(N,edo,HPlist,pulledindices,dirs,newcoords,geometry)
    
Given a lattice size `N`, a matrix encoding the aminoacids initial position `edo`, an aminoacid sequence `HPlist` ,an array containing the 
pulled indices `pulledindices`, an array containing the direction in which the chain was pulled `dirs`, a matrix containing the new positions 
for the pulled indices `newcoords`, a geometry; returns a multidimensional array containing the states at each point of the Metropolis scheme.
"""
function reconstructStates(N::Int,edo,HPlist::Vector{Amin},pulledindices,dirs::Vector{directions},newcoords,geometry::geometries)

    if geometry == fcc || geometry == triangular2D

        if geometry == fcc
            dim = 3
        else
            dim = 2
        end
        reconstructedSates = zeros(Int16,(length(HPlist),dim,length(pulledindices)+1))
        reconstructedSates[:,:,1] = edo
        for l in 2:length(pulledindices)+1
            ind = pulledindices[l-1]
    
            if ind == 0 # Pull move did not happen.
                reconstructedSates[:,:,l] = copy(reconstructedSates[:,:,l-1]) # State is unchanged
    
            else
                newedo = copy(reconstructedSates[:,:,l-1])
                refedo = copy(reconstructedSates[:,:,l-1])
                newedo[ind,:] = newcoords[l-1,:]
    
                if ind == 1
                    stateconf = validConf(N,ind,newedo,HPlist,forwards,geometry)
                    if stateconf == false
                        for k in (ind+1):length(HPlist)
                            if stateconf == false
                                newedo[k,:] = refedo[k-1,:] # New position. 
                                stateconf = validConf(N,k,newedo,HPlist,forwards,geometry) # Check whether the new configuration is valid.
                            else
                                break
                            end
                        end
                    end
                    reconstructedSates[:,:,l] = newedo
    
    
    
                elseif ind == length(HPlist)
                    stateconf = validConf(N,ind,newedo,HPlist,backwards,geometry)
                    if stateconf == false
                        for k in (ind-1):-1:1
                            if stateconf == false
                                newedo[k,:] = refedo[k+1,:] # New position. 
                                stateconf = validConf(N,k,newedo,HPlist,backwards,geometry) # Check whether the new configuration is valid.
                            else
                                break
                            end
                        end
                    end
                    reconstructedSates[:,:,l] = newedo
    
    
    
                else
                    dir = dirs[l-1]
    
                    if dir == backwards
                        stateconf = validConf(N,ind,newedo,HPlist,dir,geometry)
                        if stateconf == false
                            for k in (ind-1):-1:1
                                if stateconf == false
                                    newedo[k,:] = refedo[k+1,:] # New position. 
                                    stateconf = validConf(N,k,newedo,HPlist,dir,geometry) # Check whether the new configuration is valid.
                                else
                                    break
                                end
                            end
                        end
                        reconstructedSates[:,:,l] = newedo
    
                    elseif dir == forwards
                        stateconf = validConf(N,ind,newedo,HPlist,dir,geometry)
                        if stateconf ==  false
                            for k in (ind+1):length(HPlist)
                                if stateconf == false
                                    newedo[k,:] = refedo[k-1,:] # New position. 
                                    stateconf = validConf(N,k,newedo,HPlist,dir,geometry) # Check whether the new configuration is valid.
                                else
                                    break
                                end
                            end
                        end
                        reconstructedSates[:,:,l] = newedo
                    end
    
                    
                end
    
            end
    
        end

        return reconstructedSates


    elseif geometry == square2D
        reconstructedSates = zeros(Int16,(length(HPlist),2,length(pulledindices)+1))
        reconstructedSates[:,:,1] = edo
        for l in 2:length(pulledindices)+1
            ind = pulledindices[l-1]
    
            if ind == 0 # Pull move did not happen.
                reconstructedSates[:,:,l] = copy(reconstructedSates[:,:,l-1]) # State is unchanged
    
            else
                newedo = copy(reconstructedSates[:,:,l-1])
                refedo = copy(reconstructedSates[:,:,l-1])
                newedo[ind,:] = newcoords[1][l-1,:]
    
    
                if ind == 1
                    newedo[ind+1,:] = newcoords[2][l-1,:]
                    stateconf = validConf(N,ind+1,newedo,HPlist,forwards,geometry)
                    if stateconf == false
                        for k in (ind+2):length(HPlist)
                            if stateconf == false
                                newedo[k,:] = refedo[k-2,:] # New position. 
                                stateconf = validConf(N,k,newedo,HPlist,forwards,geometry) # Check whether the new configuration is valid.
                            else
                                break
                            end
                        end
                    end
                    reconstructedSates[:,:,l] = newedo
    
    
    
                elseif ind == length(HPlist)
                    newedo[ind-1,:] = newcoords[2][l-1,:]
                    stateconf = validConf(N,ind-1,newedo,HPlist,backwards,geometry)
                    if stateconf == false
                        for k in (ind-2):-1:1
                            if stateconf == false
                                newedo[k,:] = refedo[k+2,:] # New position.
                                stateconf = validConf(N,k,newedo,HPlist,backwards,geometry) # Check whether the new configuration is valid.
                            else
                                break
                            end
                        end
                    end
                    reconstructedSates[:,:,l] = newedo
    
    
    
                else
                    dir = dirs[l-1]
                    if dir == backwards
                        newedo[ind-1,:] = newcoords[2][l-1,:]
                        stateconf = validConf(N,ind-1,newedo,HPlist,dir,geometry)
                        if stateconf == false
                            for k in (ind-2):-1:1
                                if stateconf == false
                                    newedo[k,:] = refedo[k+2,:] # New position.
                                    stateconf = validConf(N,k,newedo,HPlist,dir,geometry) # Check whether the new configuration is valid.
                                else
                                    break
                                end
                            end
                        end
                        reconstructedSates[:,:,l] = newedo
    
                    elseif dir == forwards
                        newedo[ind+1,:] = newcoords[2][l-1,:]
                        stateconf = validConf(N,ind+1,newedo,HPlist,dir,geometry)
                        if stateconf ==  false
                            for k in (ind+2):length(HPlist)
                                if stateconf == false
                                    newedo[k,:] = refedo[k-2,:] # New position. 
                                    stateconf = validConf(N,k,newedo,HPlist,dir,geometry) # Check whether the new configuration is valid.
                                else
                                    break
                                end
                            end
                        end
                        reconstructedSates[:,:,l] = newedo
                    end
    
                end
    
            end
    
        end
    
        return reconstructedSates

    end
end
