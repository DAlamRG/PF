include("/PF/HP3D/HP-model3D-pullmoves.jl")

using Plots
    



"""
   fccLattice(n)  
Given the size of a 3D array `n`; returns a matrix whose rows are the positions on the fcc lattice..
"""
function fccLattice(n)
    a1=[1,1,0]
    a2=[0,1,1] 
    a3=[1,0,1] # Basis vectors for the Bravais lattice (fcc)
    
    fccM=zeros(Float64,(Int(n^3),3)) # This matrix will contain the linear combinations. This is the fcc lattice. 
    # Each entry fccM[i,j,k] represents the site \vec{r}=(i,k,k).
    
    for k in 1:n, j in 1:n, i in 1:n
        v=i*a1+j*a2+k*a3
        l=j+(k-1)*n
        s=i+(l-1)*n
        fccM[s,:]=v
    end
    return fccM
end














"""
   fccPositions(fccM,N,edo,HPlist)  
Given a fcc lattice `fccM` ,a 3D array size `N`, a configuration `edo`, and a sequence of aminoacids `HPlist`; returns the configuration in the fcc lattice.
"""
function fccPositions(fccM,N,edo,HPlist)
    
    edofcc=ones(Float64,(length(HPlist),3)) # This matrix will contain the positions in the fcc lattice.

    for w in 1:length(HPlist) # Iterate over the monomers.
        inds=edo[w,:]
        i,j,k=inds # Position in the previous array.
        i=mod1(i,N)
        j=mod1(j,N)
        k=mod1(k,N) # Make the indices periodic.
        l=j+(k-1)*N
        s=i+(l-1)*N # s is the fast index of the matrix `fccLattice` which returns the fcc position coordinates (`fccLattice[s,:]`).
        edofcc[w,:]=fccM[s,:]
    end

    return edofcc
end















    
    
    
# Write the visualization function. 
"""
    visHP3D(edo,HPlist,N,geometry)
Given a matrix encoding the aminoacids positions `edo`, the aminoacid sequence `HPlist`, and three sets of coordinate limits
; returns a plot displaying the protein configuration.
"""
function visHP3D(edo,HPlist,N,geometry)
    colours=[] # Assign colours to the H and P aminoacids.
    for i in 1:length(HPlist)
        el=HPlist[i]
        if el == 1
            push!(colours,"green")
        else
            push!(colours,"red")
        end
    end


    if geometry == cubic
        println("Not yet finished writing this function")

    elseif geometry == fcc
        fccM=fccLattice(N)
        edofcc=fccPositions(fccM,N,edo,HPlist)

        plt=scatter(fccM[:,1],fccM[:,2],fccM[:,3],color="gray",alpha=0.2,label="",markersize=7,grid=:false,xlabel="",
        ylabel="",zlabel="",title="Current protein configuration ")
        plot!(edofcc[:,1],edofcc[:,2],edofcc[:,3],lw=2,markershape=:circle,markercolor=colours ,markersize=7,color="purple",label="")
        xlims!((minimum(edofcc[:,1])-1,maximum(edofcc[:,1])+1))
        ylims!((minimum(edofcc[:,2])-1,maximum(edofcc[:,2])+1))
        zlims!((minimum(edofcc[:,3])-1,maximum(edofcc[:,3])+1))
        display(plt)
    end
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

