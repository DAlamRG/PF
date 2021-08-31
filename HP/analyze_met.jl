using Base: String
# This script takes the data stores in output and returns the thermodynamic analysis.

using Plots
using DelimitedFiles


include("./HP-pullmoves.jl")








"""
    analyze_met(name)
Given the name of the directory `name` where the relevant data is stored; returns the reconstructed states.
"""
function gifFolding(name::String)
 
    pathname1="./output"*"/"*name*"/"
    # pathname2="./output"*"/"*name*"/"*string(nrun)*"_"
    
    temperatures=readdlm(pathname1*"temperatures.csv",',')
    HPlist=readdlm(pathname1*"HPlist.csv",',')
    initialconf=readdlm(pathname1*"initialconf.csv",',')
    N=Int(readdlm(pathname1*"latticesize.csv",',')[1])
    nums=Int(readdlm(pathname1*"mc_sweeps.csv",',')[1])
    geometry=geom_int(readdlm(pathname1*"geometry.csv",',')[1])
    pfmodel=readdlm(pathname1*"pfmodel.csv",',')


    ns=nums*length(HPlist) # Number of total iterations in the simulation.
    ft=Int((length(temperatures)*ns)+1) # Need to add the second term beacuse I am storing the initial configuration.
   
    # Now I need to reconstruct the visited states, according to the geometry
    dim=length(initialconf[1,:])

    edos=zeros(Int64,(length(HPlist),dim,ft))
    pulledindices=Int.(readdlm(pathname2*"1_1.csv",','))
    dirs=readdlm(pathname2*"1_2.csv",',')
    dirs=dirsf(dirs) # Data is written as Int type, turn it back to enum type.

    if geometry == triangular2D || geometry == fcc
        newcoords=readdlm(pathname2*"1_3.csv",',')
        rd=reconstructStates(N,initialconf,HPlist,pulledindices,dirs,newcoords,geometry)
        edos[:,:,1:ns+1]=rd
        laststate=rd[:,:,end]
        
        cont=ns+2
        for i in 2:length(temperatures)
            pulledindices=Int.(readdlm(pathname2*string(i)*"_1.csv",','))
            dirs=readdlm(pathname2*string(i)*"_2.csv",',')
            dirs=dirsf(dirs)
            newcoords=readdlm(pathname2*string(i)*"_3.csv",',')
            rS=reconstructStates(N,laststate,HPlist,pulledindices,dirs,newcoords,geometry)
            edos[:,:,cont:cont+(ns-1)]=rS[:,:,2:end] 
            laststate=rS[:,:,end]
            cont=cont+ns
        end


    elseif geometry == square2D
        newcoords1=readdlm(pathname2*"1_3_1.csv",',')
        newcoords2=readdlm(pathname2*"1_3_2.csv",',')
        newcoords=(newcoords1,newcoords2)
        rd=reconstructStates(N,initialconf,HPlist,pulledindices,dirs,newcoords,geometry)
        edos[:,:,1:ns+1]=rd
        laststate=rd[:,:,end]

        cont=ns+2
        for i in 2:length(temperatures)
            pulledindices=Int.(readdlm(pathname2*string(i)*"_1.csv",','))
            dirs=readdlm(pathname2*string(i)*"_2.csv",',')
            dirs=dirsf(dirs)
            newcoords1=readdlm(pathname2*string(i)*"_3_1.csv",',')
            newcoords2=readdlm(pathname2*string(i)*"_3_2.csv",',')
            newcoords=(newcoords1,newcoords2)
            rS=reconstructStates(N,laststate,HPlist,pulledindices,dirs,newcoords,geometry)
            edos[:,:,cont:cont+(ns-1)]=rS[:,:,2:end] 
            laststate=rS[:,:,end]
            cont=cont+ns
        end
    end

   
end
