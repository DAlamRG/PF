using Plots: length, get_fillalpha
using Base: Float64, Int16, String

# This script takes the data stores in output and returns the thermodynamic analysis.

using Plots
gr()
using DelimitedFiles
using Statistics
using LinearAlgebra


include("./HP-pullmoves.jl")








"""
    center_mass(vec,dim)
Given a matrix encoding the monomers' positions `vec` and a dimension; returns a the center of mass vector.
"""
function center_mass(vec,dim)
    m = size(vec)[1]
    rcm = zeros(Float64,dim)
    for d in 1:dim
        rcm[d] = sum(vec[:,d])/m
    end
    return rcm
end









"""
    Rᵧ(vec)
Given an array encoding the positions for the monomers on the chain `vec`; returns the corresponding radius of giration. 
"""
function Rᵧ(vec)
    dim = size(vec)[2] 
    m = size(vec)[1]
    rcm = center_mass(vec,dim)
    Rg = 0.0
    for i in 1:m
        vaux = vec[i,:]-rcm
        Rg = Rg + norm(vaux)^2
    end
    Rg = sqrt(Rg/m)
    return Rg
end










"""
    ee_distance(vec)
Given a mtric encoding the aminoacid's positions; returns the distance between the first and last monomers.
"""
function ee_distance(vec)
    ree = norm(vec[1,:]-vec[end,:])
    return ree
end










"""
    states_met_run(name,nrun)
Given the name of the directory `name` where the relevant data is stored and the number of the independent run `nrun` to be reconstructed; returns the 
reconstructed states for the given run.
"""
function states_met_run(name::String,nrun::Int)

    # Load the data.
    pathname1 = "./output"*"/"*name*"/"

    geometry = geom_int(Int(readdlm(pathname1*"geometry.csv",',')[1]))
    HPlist = readdlm(pathname1*"HPlist.csv",',')
    HPlist = Amin[amin_dict[j] for j in HPlist]
    initialconf = Int16.(readdlm(pathname1*"initialconf.csv",','))
    N = Int(readdlm(pathname1*"latticesize.csv",',')[1])
    nums = readdlm(pathname1*"mc_sweeps.csv",',')[1]
    temperatures = readdlm(pathname1*"temperatures.csv",',')


    n_amin = length(HPlist) 
    ns = Int(nums*n_amin)
    ft = Int(length(temperatures)*ns) # Total number of visted states, minus the initial one.

    # Now I need to reconstruct the visited states, according to the geometry
    dim = length(initialconf[1,:])
    states = zeros(Int16,(n_amin,dim,ft))
    
    pathname2 = "./output"*"/"*name*"/"*string(nrun)*"_"
    pulledindices = Int.(readdlm(pathname2*"1_1.csv",','))
    dirs = dirsf(readdlm(pathname2*"1_2.csv",',')) # Data is written as Int type, turn it back to enum type.

    if geometry == triangular2D || geometry == fcc
        newcoords = readdlm(pathname2*"1_3.csv",',')
        rd = reconstructStates(N,initialconf,HPlist,pulledindices,dirs,newcoords,geometry)
        states[:,:,1:ns] = rd[:,:,2:end]
        laststate = rd[:,:,end]
        
        cont = ns+1
        for i in 2:length(temperatures)
            pulledindices = Int.(readdlm(pathname2*string(i)*"_1.csv",','))
            dirs = dirsf(readdlm(pathname2*string(i)*"_2.csv",','))
            newcoords = readdlm(pathname2*string(i)*"_3.csv",',')
            rd = reconstructStates(N,laststate,HPlist,pulledindices,dirs,newcoords,geometry)
            states[:,:,cont:cont+(ns-1)] = rd[:,:,2:end] 
            laststate = rd[:,:,end]
            cont = cont + ns
        end


    elseif geometry == square2D
        newcoords1 = readdlm(pathname2*"1_3_1.csv",',')
        newcoords2 = readdlm(pathname2*"1_3_2.csv",',')
        newcoords = (newcoords1,newcoords2)
        rd = reconstructStates(N,initialconf,HPlist,pulledindices,dirs,newcoords,geometry)
        states[:,:,1:ns] = rd[:,:,2:end]
        laststate = rd[:,:,end]

        cont = ns+1
        for i in 2:length(temperatures)
            pulledindices = Int.(readdlm(pathname2*string(i)*"_1.csv",','))
            dirs = readdlm(pathname2*string(i)*"_2.csv",',')
            dirs = dirsf(dirs)
            newcoords1 = readdlm(pathname2*string(i)*"_3_1.csv",',')
            newcoords2 = readdlm(pathname2*string(i)*"_3_2.csv",',')
            newcoords = (newcoords1,newcoords2)
            rd = reconstructStates(N,laststate,HPlist,pulledindices,dirs,newcoords,geometry)
            edos[:,:,cont:cont+(ns-1)] = rd[:,:,2:end] 
            laststate = rd[:,:,end]
            cont = cont+ns
        end
    end


    return states
end











"""
    states_met(name)
Given a directory name `name`; returns an array containing the visited states during all of the runs.
"""
function states_met(name::String)

    # Load the data.
    pathname1 = "./output"*"/"*name*"/"

    nruns = Int(readdlm(pathname1*"nruns.csv",',')[1])
    geometry = geom_int(Int(readdlm(pathname1*"geometry.csv",',')[1]))
    HPlist = readdlm(pathname1*"HPlist.csv",',')
    HPlist = Amin[amin_dict[j] for j in HPlist]
    initialconf = Int16.(readdlm(pathname1*"initialconf.csv",','))
    N = Int(readdlm(pathname1*"latticesize.csv",',')[1])
    nums = readdlm(pathname1*"mc_sweeps.csv",',')[1]
    pfmodel = readdlm(pathname1*"pfmodel.csv",',')
    temperatures = readdlm(pathname1*"temperatures.csv",',')

    n_amin = length(HPlist) 
    ns = Int(nums*n_amin)
    ft = Int(length(temperatures)*ns) # Total number of visted states, minus the initial one.
    dim = length(initialconf[1,:])


    run_states = fill(zeros(Int16,(n_amin,dim,ft)),nruns) # This vector will contain all of the states for each run. We need to fill it.

    for l in 1:nruns
        run_states[l] = states_met_run(name,l)
    end

    return run_states

end















"""
    analyze_met_thermo(name)
Given the name of the directory `name` where the relevant data is stored; stores the values for the internal energies and heat capacity, 
to the same directory.
"""
function analyze_met_thermo(name::String)

    # Load the data.
    pathname1 = "./output"*"/"*name*"/"

    nruns = Int(readdlm(pathname1*"nruns.csv",',')[1])
    geometry = geom_int(Int(readdlm(pathname1*"geometry.csv",',')[1]))
    HPlist = readdlm(pathname1*"HPlist.csv",',')
    HPlist = Amin[amin_dict[j] for j in HPlist]
    initialconf = readdlm(pathname1*"initialconf.csv",',')
    N = Int(readdlm(pathname1*"latticesize.csv",',')[1])
    nums = readdlm(pathname1*"mc_sweeps.csv",',')[1]
    pfmodel = readdlm(pathname1*"pfmodel.csv",',')
    temperatures = readdlm(pathname1*"temperatures.csv",',')



    # Now, we need to retrieve the energies for all the visited temperatures, for all the independent runs.
    n_amin = length(HPlist) 
    ns = Int(nums*n_amin)
    n_middle = Int(ceil(ns/2))

    us = zeros(Float64,(length(temperatures),nruns)) # This matrix will contain "half" the energies at each temperature, and for each independent run.
    cs = zeros(Float64,(length(temperatures),nruns))

    for l in 1:nruns
        # First, import the energies.
        pathnameaux = pathname1*"/"*string(l)*"_energies.csv"
        run_energies = readdlm(pathnameaux,',')[2:end]

        cont = 1
        for k in length(temperatures):-1:1
            es = run_energies[cont+n_middle:cont+(ns-1)]
            uk = mean(es)/n_amin
            us[k,l] = uk

            β = 1/temperatures[k]
            cs[k,l] = ((β^2)/n_amin)*(mean(es.^2)-(mean(es)^2))
            cont = cont + ns
        end
        println("Finished with run number: ",l)
    end

    # Save the thermodynamic variables.
    writedlm(pathname1*"us.csv",us,',')
    writedlm(pathname1*"cs.csv",cs,',')
end

# analyze_met_thermo("simu3")



################   This bit is to compare the energies between the square and triangular geometries    ################
# pathnameEnergy1 = "./output"*"/"*"simu4"*"/1_energies.csv"

# pathnameEnergy2 = "./output"*"/"*"simu5"*"/1_energies.csv"

# ens1 = readdlm(pathnameEnergy1,',')
# ens2 = readdlm(pathnameEnergy2,',')

# plot(1:length(ens1),ens1,color="green",lw=2,xlabel="step",ylabel="Energy per monomer",label="",alpha=0.8)
# plot!(1:length(ens2),ens2,lw=1.5,alpha=0.5,label="",color="purple")





############# This bit is here to analyze the third simulation ###############

us3 = readdlm("./output/simu3/us.csv",',')
cs3 = readdlm("./output/simu3/cs.csv",',')
temps = readdlm("./output/simu3/temperatures.csv",',')
up = Float64[mean(us3[k,:]) for k in 1:length(temps)]
uσs = Float64[std(us3[k,:]) for k in 1:length(temps)]
cp = Float64[mean(cs3[k,:]) for k in 1:length(temps)]
cσs = Float64[std(cs3[k,:]) for k in 1:length(temps)]


# display(scatter(temps,us3[:,1],ms=2,color="green",alpha=0.7,xlabel="T",ylabel="u(T)",title="Internal energy per monomer for simu3",label=""))
# display(scatter(temps,cs3[:,1],ms=3,color="purple",alpha=0.7,xlabel="T",ylabel="c(T)",title="heat per monomer for simu3",label=""))

display(scatter(temps,up,ms=3.5,color="green",alpha=0.7,xlabel="T",ylabel="u(T)",title="Internal energy per monomer for simu3",label="",ribbon=uσs,fillalpha=0.3))
# display(scatter(temps,cp,ms=3.5,color="purple",alpha=0.7,xlabel="T",ylabel="c(T)",title="Specific heat per monomer for simu3",label="",ribbon=cσs,fillalpha=0.3))