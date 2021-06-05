# I write a function which returns a plot which allows for better visualization of the amino acids type, positions and bonds.
using Plots






"""
    triangularLattice(nα,nβ)
Returns a matrix with the coordinates for `nα*nβ` points in a triangular lattice.
"""
function triangularLattice(nα,nβ)
    v1=[1,0]
    v2=[-1/2,sqrt(3)/2] # Basis vectors for the triangular Bravais lattice.
    
    bM=zeros(Float64,(Int(nα*nβ),2)) # This matrix will contain the linear combinations.
    
    for β in 1:nβ, α in 1:nα
        v=α*v1+β*v2
        k=α+(β-1)*nα
        bM[k,:]=v
    end
    return bM
end











"""
    chBasisTriangular(m)
Given a matrix `m` whose rows are the coordinates in a regular square lattice; returns the same elements but in the triangular lattice.
"""
function chBasisTriangular(m)
    A=[[-1/2 1];[sqrt(3)/2 0]] # CHange of basis matrix.
    mp=zeros(Float64,size(m)) # This array contains the vectors expressed in the new basis.
    for k in 1:size(m)[1]
        el=m[k,:]
        nv=A*el
        mp[k,:]=nv
    end
    return mp
end












# Write the visualization function. This one takes the geometry into account.
"""
    visHP2D(edo,HPlist,nα,nβ,geometry)
Given a matrix encoding the aminoacids positions `edo`, the aminoacid sequence `HPlist`, two numbers `nα,nβ` representing the limits in 
the `x` and `y` axis, and a geometry; returns a plot displaying the protein configuration.
"""
function visHP2D(edo,HPlist,nα,nβ,geometry)
    colours=[] # Contains the colors for the two kinds of aminoacids (H = -1 = "red", P = 1 = "green").
    for i in 1:length(HPlist) 
        el=HPlist[i]
        if el == 1
            push!(colours,"green")
        else
            push!(colours,"red")
        end
    end



    if geometry == square2D
        
        plt=plot(edo[:,2],edo[:,1],lw=2,markershape=:circle,markercolor=colours ,markersize=7,color="purple",label="",xlabel="",
        ylabel="",title="Current protein configuration ")
        xlims!((0,nα))
        ylims!((0,nβ))

    elseif geometry == triangular2D
        tl=triangularLattice(nα,nβ) # Generate the background lattice
        edo=chBasisTriangular(edo) # Change the basis for the current configuration `edo`.

        plt= scatter(tl[:,1],tl[:,2],color="gray",alpha=0.2,label="",markersize=7,grid=:false,xlabel="",
        ylabel="",title="Current protein configuration ")
        plot!(edo[:,1],edo[:,2],lw=2,markershape=:circle,markercolor=colours ,markersize=7,color="purple",label="")
        xlims!((-(0.5)*nβ,nα))
        ylims!((0,nβ*0.866025))
    end
    display(plt)
end
















# 
"""
    visHP2Dheat(red,edo)

Given the array `red` containing the protein´s configuration, a matrix encoding the aminoacids positions `edo`; returns 
a heat map with the protein´s configuration.
"""
function visHP2Dheat(red,edo)
    N=size(red)[1]
    plt=heatmap(1:N,1:N,red,xlabel="",ylabel="",label="",title="Current protein configuration")
    plot!(edo[:,2],edo[:,1],lw=2,markershape=:circle,markersize=3,color="purple",label="")
    display(plt)
end












# I would also like to be able to visulaize the protein´s energies across iterations of the Metropolis scheme.
"""
        energiesHP(T,energies,HPlist)
    
Given the temperature `T`, an array containing the visited energies during the simulation `energies`, and the sequence of aminoacids
`HPlist`; returns a plot of the energies as a function of iteration step.
"""
function energiesHP(T,energies,HPlist)
    n=(length(energies)-1)/length(HPlist)
    pp=plot(range(1,length=(length(energies)-1),stop=n),energies[1:end-1],color="green",lw=2,alpha=0.7,xlabel="iteration per site",
    ylabel="E",label="",title="Visited energies at T=$T")
    display(pp)
end
















# I would also like to have the energies histogram.
"""
        energiesdistHP(T,energies)
    
Given the temperature `T`, and an array containing the visited energies during the simulation `energies`; returns 
a histogram of the visited energies (vertical axis is in logarithmic scale).
"""
function energiesdistHP(T,energies)
    n=Int(abs(minimum(energies)))
    enp=energies.-1
    pp=histogram(enp,nbins=n,yaxis=:log,color="green",title="Histogram of visited energies at T=$T",alpha=0.7,
    label="",xlabel="E",ylabel="log(n(E))")
    display(pp)
end