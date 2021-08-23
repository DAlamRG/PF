@enum lattAmin::Int8 begin
    Null = 0
    R = 1
    D = 2
    N = 3
    A = 4
end

const HPXN_Ps = [R]
const HPXN_Ns = [D]
const HPXN_Xs = [N]
const HPXN_Hs = [A]

function HPXNgroupNumber(amin::lattAmin)::Int8
    if amin in HPXN_Ps
        return 1
    elseif amin in HPXN_Ns
        return 2
    elseif amin in HPXN_Xs
        return 3
    elseif amin in HPXN_Hs
        return 4
    else
        return 5
    end
end

const HPXNinterMatrix = [
    [-4 0 0 0 0] ; 
    [0 1 -1 0 0] ; 
    [0 -1 1 0 0] ; 
    [0 0 0 0 0] ; 
    [0 0 0 0 0] 
]

function interaction(amin1::lattAmin,amin2::lattAmin,groupFunction::Function,int_matrix::Matrix{<:Real})::Float64
    i = groupFunction(amin1)
    j = groupFunction(amin2)
    return int_matrix[i,j]
end
