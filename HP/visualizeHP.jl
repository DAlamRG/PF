using Base: String, Int16, Int8
using Plots: length, push!
# This script contains the code necessary for visualizing the outputs of the main function.


using Plots
using DelimitedFiles
using Colors
pyplot()



include("./HP-pullmoves.jl")











"""
    triangularLattice(nα,nβ)
Returns a matrix with the coordinates for `nα*nβ` points in a triangular lattice.
"""
function triangularLattice(nα,nβ)
    v1 = [1,0]
    v2 = [-1/2,sqrt(3)/2] # Basis vectors for the triangular Bravais lattice.
    
    bM = zeros(Float64,(Int(nα*nβ),2)) # This matrix will contain the linear combinations.
    
    for β in 1:nβ, α in 1:nα
        v = α*v1 + β*v2
        k = α + (β-1)*nα
        bM[k,:] = v
    end
    return bM
end















"""
    chBasisTriangular(m)
Given a matrix `m` whose rows are the coordinates in a regular square lattice; returns the same elements but in the triangular lattice.
"""
function chBasisTriangular(m)
    A = [[-1/2 1];[sqrt(3)/2 0]] # Change of basis matrix.
    mp = zeros(Float64,size(m)) # This array contains the vectors expressed in the new basis.
    for k in 1:size(m)[1]
        el = m[k,:]
        nv = A*el
        mp[k,:] = nv
    end
    return mp
end














"""
   fccLattice(n)  
Given the size of a 3D array `n`; returns a matrix whose rows are the positions on the fcc lattice..
"""
function fccLattice(n)
    a1 = [1,1,0]
    a2 = [0,1,1] 
    a3 = [1,0,1] # Basis vectors for the Bravais lattice (fcc)
    
    fccM = zeros(Float64,(Int(n^3),3)) # This matrix will contain the linear combinations. This is the fcc lattice. 
    # Each entry fccM[i,j,k] represents the site \vec{r}=(i,k,k).
    
    for k in 1:n, j in 1:n, i in 1:n
        v = i*a1 + j*a2 + k*a3
        l = j + (k-1)*n
        s = i + (l-1)*n
        fccM[s,:] = v
    end
    return fccM
end














"""
   fccPositions(fccM,N,edo,HPlist)  
Given a fcc lattice `fccM` ,a 3D array size `N`, a configuration `edo`, and a sequence of aminoacids `HPlist`; returns the configuration in the fcc lattice.
"""
function fccPositions(fccM,N,edo,HPlist)
    
    edofcc = ones(Float64,(length(HPlist),3)) # This matrix will contain the positions in the fcc lattice.

    for w in 1:length(HPlist) # Iterate over the monomers.
        inds = edo[w,:]
        i,j,k = inds # Position in the previous array.
        i = mod1(i,N)
        j = mod1(j,N)
        k = mod1(k,N) # Make the indices periodic.
        l = j + (k-1)*N
        s = i + (l-1)*N # s is the fast index of the matrix `fccLattice` which returns the fcc position coordinates (`fccLattice[s,:]`).
        edofcc[w,:] = fccM[s,:]
    end

    return edofcc
end












"""
    color_amin(HPlist)

Given a list of aminoacids `HPlist`, returns a list ofthe corresponding colours. 
"""
function color_amin(HPlist)
    colours = String[]

    for i in 1:length(HPlist)
        el = HPlist[i]
        if el == h
            push!(colours,"crimson")
        elseif el == H
            push!(colours,"red")
        elseif el == P
            push!(colours,"green")
        elseif el == N
            push!(colours,"lightgreen")
        elseif el == X
            push!(colours,"turquoise")
        elseif el == Y
            push!(colours,"darkorange1")
        end
    end
    return colours
end












"""
    visHP(edo,HPlist,N,geometry,time,temp,camargs)
Given a matrix encoding the aminoacids positions `edo`, the aminoacid sequence `HPlist`, three sets of coordinate limits, a 
time `time`, a temperature `temp`, and the arguments for the camera `camargs`; returns a plot displaying the protein configuration.
"""
function visHP(edo,HPlist,N,geometry,time,temp,camargs)
    colours = color_amin(HPlist)
    

    if geometry == square2D
        
        plt = plot(edo[:,2],edo[:,1],lw=2,markershape=:circle,markercolor=colours ,markersize=7,color="purple",label="",xlabel="",
        ylabel="",title="Current protein configuration \n (Time,Temperature)=($time , $temp) ")
        xlims!((0,N))
        ylims!((0,N))

    elseif geometry == triangular2D
        tl = triangularLattice(N,N) # Generate the background lattice
        edo = chBasisTriangular(edo) # Change the basis for the current configuration `edo`.

        plt= scatter(tl[:,1],tl[:,2],color="gray",alpha=0.2,label="",markersize=7,grid=:false,xlabel="",
        ylabel ="",title="Current protein configuration \n (Time,Temperature)=($time , $temp) ")
        plot!(edo[:,1],edo[:,2],lw=2,markershape=:circle,markercolor=colours ,markersize=7,color="purple",label="")
        xlims!((-(0.5)*N,N))
        ylims!((0,N*0.866025))
    
    elseif geometry == cubic
        println("Work in progress")
    elseif geometry == fcc
        fccM = fccLattice(N)
        edofcc = fccPositions(fccM,N,edo,HPlist)

        if camargs[1] 
            plt = scatter(fccM[:,1],fccM[:,2],fccM[:,3],color="gray",alpha=0.01,label="",markersize=7,camera=camargs[2]
            ,grid=:false,xlabel="",ylabel="",zlabel=""
            ,title="Current protein configuration \n (Time,Temperature)=($time , $temp)")
            plot!(edofcc[:,1],edofcc[:,2],edofcc[:,3],lw=2,markershape=:circle,markercolor=colours ,camera=camargs[2],
            markersize=7,color="purple",label="")
            xlims!((minimum(edofcc[:,1])-1,maximum(edofcc[:,1])+1))
            ylims!((minimum(edofcc[:,2])-1,maximum(edofcc[:,2])+1))
            zlims!((minimum(edofcc[:,3])-1,maximum(edofcc[:,3])+1))


        else
            plt = scatter(fccM[:,1],fccM[:,2],fccM[:,3],color="gray",alpha=0.2,label="",markersize=7,grid=:false,xlabel="",
            ylabel="",zlabel="",title="Current protein configuration \n (Time,Temperature)=($time , $temp)")
            plot!(edofcc[:,1],edofcc[:,2],edofcc[:,3],lw=2,markershape=:circle,markercolor=colours ,markersize=7,color="purple",label="")
            xlims!((minimum(edofcc[:,1])-1,maximum(edofcc[:,1])+1))
            ylims!((minimum(edofcc[:,2])-1,maximum(edofcc[:,2])+1))
            zlims!((minimum(edofcc[:,3])-1,maximum(edofcc[:,3])+1))

        end

    end

    return plt

end
        
















"""
    gifFolding(name,nrun,nskip,gifname)
Given he name of the directory where the data is stored `name`, the run number in the data `nrun`,  the number of frames to be 
skipped in the animation `nskip`, and a name for the animation `gifname`; saves an animation of the folding process to the directory from where 
the data comes from.
"""
function gifFolding(name,nrun,nskip,gifname)
 
    pathname1 = "./output"*"/"*name*"/"
    pathname2 = "./output"*"/"*name*"/"*string(nrun)*"_"
    
    temperatures = readdlm(pathname1*"temperatures.csv",',')
    HPlist = readdlm(pathname1*"HPlist.csv",',')
    initialconf = readdlm(pathname1*"initialconf.csv",',')
    N = Int(readdlm(pathname1*"latticesize.csv",',')[1])
    nums = Int(readdlm(pathname1*"mc_sweeps.csv",',')[1])
    geometry = geom_int(readdlm(pathname1*"geometry.csv",',')[1])


    ns = nums*length(HPlist) # Number of total iterations in the simulation.
    ft = Int((length(temperatures)*ns)+1) # Need to add the second term beacuse I am storing the initial configuration.
    times = 1:nskip:ft # Number of plots to be made.

    # Now I need to reconstruct the visited states, according to the geometry
    dim = length(initialconf[1,:])

    edos = zeros(Int64,(length(HPlist),dim,ft))
    pulledindices = Int.(readdlm(pathname2*"1_1.csv",','))
    dirs = readdlm(pathname2*"1_2.csv",',')
    dirs = dirsf(dirs) # Data is written as Int type, turn it back to enum type.

    if geometry == triangular2D || geometry == fcc
        newcoords = readdlm(pathname2*"1_3.csv",',')
        rd = reconstructStates(N,initialconf,HPlist,pulledindices,dirs,newcoords,geometry)
        edos[:,:,1:ns+1] = rd
        laststate = rd[:,:,end]
        
        cont = ns+2
        for i in 2:length(temperatures)
            pulledindices = Int.(readdlm(pathname2*string(i)*"_1.csv",','))
            dirs = readdlm(pathname2*string(i)*"_2.csv",',')
            dirs = dirsf(dirs)
            newcoords = readdlm(pathname2*string(i)*"_3.csv",',')
            rS = reconstructStates(N,laststate,HPlist,pulledindices,dirs,newcoords,geometry)
            edos[:,:,cont:cont+(ns-1)] = rS[:,:,2:end] 
            laststate = rS[:,:,end]
            cont = cont+ns
        end


    elseif geometry == square2D
        newcoords1 = readdlm(pathname2*"1_3_1.csv",',')
        newcoords2 = readdlm(pathname2*"1_3_2.csv",',')
        newcoords = (newcoords1,newcoords2)
        rd = reconstructStates(N,initialconf,HPlist,pulledindices,dirs,newcoords,geometry)
        edos[:,:,1:ns+1] = rd
        laststate = rd[:,:,end]

        cont = ns+2
        for i in 2:length(temperatures)
            pulledindices = Int.(readdlm(pathname2*string(i)*"_1.csv",','))
            dirs = readdlm(pathname2*string(i)*"_2.csv",',')
            dirs = dirsf(dirs)
            newcoords1 = readdlm(pathname2*string(i)*"_3_1.csv",',')
            newcoords2 = readdlm(pathname2*string(i)*"_3_2.csv",',')
            newcoords = (newcoords1,newcoords2)
            rS = reconstructStates(N,laststate,HPlist,pulledindices,dirs,newcoords,geometry)
            edos[:,:,cont:cont+(ns-1)] = rS[:,:,2:end] 
            laststate = rS[:,:,end]
            cont = cont+ns
        end
    end

    println("States have been reconstructed.")

    temps = reverse(temperatures)
    if dim == 2
        animfold = @animate for i in times 
            mcnumber = Int(ceil(((i-1)/length(HPlist)))) # Monte Carlo step.
            ind = 0
            if i == 1
                ind = 1
            else
                ind = Int(ceil((i-1)/ns))
            end
            temp = temps[ind]
            visHP(edos[:,:,i],HPlist,N,geometry,mcnumber,round(temp,digits=2),(false,(0,0)))
        end

    elseif dim == 3
        counteranim = 1
        animfold = @animate for i in times
            counteranim += 1 
            mcnumber = Int(ceil(((i-1)/length(HPlist)))) # Monte Carlo step.
            ind = 0
            if i == 1
                ind = 1
            else
                ind = Int(ceil((i-1)/ns))
            end
            temp = temps[ind]
            cameraval = (mod1(2*counteranim,360),30)
            visHP(edos[:,:,i],HPlist,N,geometry,mcnumber,round(temp,digits=2),(true,cameraval))
        end
    end
    println("Animation is stored, all that is left to do is to save it to a .gif file.")
    gif(animfold,pathname1*gifname,fps=24)

end











# Test (Ahhh!)

gifFolding("simu1",1,3,"simu2_0.gif")
