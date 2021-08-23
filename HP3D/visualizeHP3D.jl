include("./HP-model3D-pullmoves.jl")

using Plots
using DelimitedFiles

pyplot()
    



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
    visHP3D(edo,HPlist,N,geometry,time,temp,camargs)
Given a matrix encoding the aminoacids positions `edo`, the aminoacid sequence `HPlist`, three sets of coordinate limits, a 
time `time`, a temperature `temp`, and the arguments for the camera `camargs`; returns a plot displaying the protein configuration.
"""
function visHP3D(edo,HPlist,N,geometry,time,temp, camargs)
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

        if camargs[1] 
            plt=scatter(fccM[:,1],fccM[:,2],fccM[:,3],color="gray",alpha=0.01,label="",markersize=7,camera=camargs[2]
            ,grid=:false,xlabel="",ylabel="",zlabel=""
            ,title="Current protein configuration \n (Time,Temperature)=($time , $temp)")
            plot!(edofcc[:,1],edofcc[:,2],edofcc[:,3],lw=2,markershape=:circle,markercolor=colours ,camera=camargs[2],
            markersize=7,color="purple",label="")
            xlims!((minimum(edofcc[:,1])-1,maximum(edofcc[:,1])+1))
            ylims!((minimum(edofcc[:,2])-1,maximum(edofcc[:,2])+1))
            zlims!((minimum(edofcc[:,3])-1,maximum(edofcc[:,3])+1))
            return plt

        else
            plt=scatter(fccM[:,1],fccM[:,2],fccM[:,3],color="gray",alpha=0.2,label="",markersize=7,grid=:false,xlabel="",
            ylabel="",zlabel="",title="Current protein configuration \n (Time,Temperature)=($time , $temp)")
            plot!(edofcc[:,1],edofcc[:,2],edofcc[:,3],lw=2,markershape=:circle,markercolor=colours ,markersize=7,color="purple",label="")
            xlims!((minimum(edofcc[:,1])-1,maximum(edofcc[:,1])+1))
            ylims!((minimum(edofcc[:,2])-1,maximum(edofcc[:,2])+1))
            zlims!((minimum(edofcc[:,3])-1,maximum(edofcc[:,3])+1))
            return plt
        end
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

















"""
    dirsf(vec)
Given a vector with integer entries; returns the equivalent `dirs` vector.
"""
function dirsf(vec)
    l=length(vec)
    dirsVec=Vector{directions}(undef,l) # Declare the equivalent `dirs` vector.
    for k in 1:l
        val=vec[k]
        if val == 1
            dirsVec[k]=forwards
        else
            dirsVec[k]=backwards
        end
    end
    return dirsVec
end














# Write a script to animate the protein folding process.

"""
    gifFolding(N,name,nums,nrun,geometry,nskip,gifname)
Given a lattice size `N`, the name where the data is stored `name`, the number of monte-carlo sweeps per temperature `nums` for
the provided data, the run number in the data `nrun`, the geomtery of the data `geometry`, the number of frames to be 
skipped in the animation `nskip`, and a name for the animation `gifname`; saves an animation to the directory from where 
the data comes from.
"""
function gifFolding(N,name,nums,nrun,geometry,nskip,gifname)
    pathname1="./output3D"*"/"*name*"/"
    pathname2="./output3D"*"/"*name*"/"*string(nrun)*"_"

    temperatures=readdlm(pathname1*"temperatures.csv",',')
    HPlist=readdlm(pathname1*"HPlist.csv",',')
    initialconf=readdlm(pathname1*"initialconf.csv",',')

    ns=nums*length(HPlist)
    ft=(length(temperatures)*ns)+length(temperatures) # Need to add the second term beacuse I am storing the final/initial configuration.
    times=1:nskip:ft # Number of plots to be made.

    # Now I need to reconstruct the visited states, according to the geometry
    edos=zeros(Int64,(length(HPlist),3,ft))
    pulledindices=Int.(readdlm(pathname2*"1_1.csv",','))
    dirs=readdlm(pathname2*"1_2.csv",',')
    dirs=dirsf(dirs) # Data is written as Int type, turn it back to enum type.

    newcoords=readdlm(pathname2*"1_3.csv",',')
    rd=reconstructStates3D(N,initialconf,HPlist,pulledindices,dirs,newcoords,geometry)
    edos[:,:,1:ns+1]=rd
    laststate=rd[:,:,end]
    cont=ns+2

    for i in 2:length(temperatures)
        pulledindices=Int.(readdlm(pathname2*string(i)*"_1.csv",','))
        dirs=readdlm(pathname2*string(i)*"_2.csv",',')
        dirs=dirsf(dirs)
        newcoords=readdlm(pathname2*string(i)*"_3.csv",',')
        rS=reconstructStates3D(N,laststate,HPlist,pulledindices,dirs,newcoords,geometry)
        edos[:,:,cont:cont+ns]=rS
        laststate=rS[:,:,end]
        cont=cont+ns+1
    end


    temps=reverse(temperatures)
    counteranim = 1
    animfold = @animate for i in times
        counteranim += 1 
        mcnumber= Int(ceil((i/length(HPlist))))
        ind=Int(ceil((i/(ns+1))))
        temp=temps[ind]
        cameraval=(mod1(2*counteranim,360),40)
        visHP3D(edos[:,:,i],HPlist,N,geometry,mcnumber,round(temp,digits=2),(true,cameraval))
    end
    println("animation is stored, all that is left to do is to save it to a .gif file")
    gif(animfold,pathname1*gifname,fps=24)
end



# Example

 gifFolding(20,"trial3Dfcc3",8,1,fcc,3,"3Danim1.gif")