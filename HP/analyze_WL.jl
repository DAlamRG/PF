
using Base: Float64, Int16, Int8
# This script analyzes the data stored in outputWL

using Plots
using StatsBase
using Statistics
using DelimitedFiles
gr()


include("./Energy.jl")





"""
    determineλ(lngE,T)
Given a dictionary `lngE` whose keys(values) are energies(log of energy densities), and a temperature `T`; returns the largest value
λ of the differences between energy density `ln(g(E))` and `βE` (i.e. `max{log(g(Eᵢ))-βE}`).
"""
function determineλ(lngE::Dict{Float64,Float64},T::Float64)
    β = 1/T
    difs = Float64[]
    for el in lngE
        e,lnge = el 
        dif = lnge-(β*e)
        push!(difs,dif)
    end
    λ = maximum(difs)
    return λ
end














"""
    sExp(E,T,lnge,λ)
Given a temperature `T`, an energy `E`, the logarithm of the enegy density `lnge` and a number `λ`; returns e^{ln(g(e))-βE-λ}.
"""
function sExp(E,T,lnge,λ)
    β = 1/T
    res = exp(lnge-β*E-λ)
    return res
end












"""
    thermo_WL(ti,tf,nTs,name)
Given an initial/final temperature `ti/tf`, a number of temperatures `nTs` and the name for the directory where the data is stored
`name`; returns a tuple containing the thermodynamic variables' values at each temperature.
"""
function thermo_WL(ti,tf,nTs,name::String)
    
    temperatures = range(ti,stop=tf,length=nTs)

    # Obtain the data path name.
    pathstring="./outputWL/"
    pathname = pathstring*name
    lngE = readdlm(pathname*"/lngE.csv",',')
    HPlist = vec(readdlm(pathname*"/HPlist.csv",','))
    n = length(HPlist)

    lngE = Dict{Float64,Float64}(lngE[i,1] => lngE[i,2] for i in 1:length(lngE[:,1]))

    # Declare arrays which will contain thermodynamic data.
    us = ones(Float64,nTs) 
    cs = ones(Float64,nTs) 
    Fs = ones(Float64,nTs)
    Ss = ones(Float64,nTs)
    
    for i in 1:length(temperatures)
        T = temperatures[i]
        𝚬 = 0 # Average energy.
        𝚬sq = 0 # Average of the squared energy.
        𝚭 = 0 # Partition function, normalizes the average energy.
        λ = determineλ(lngE,T)
        for element in lngE
            e,lnge = element # Extract the energy and natural logarithm of the energy density.
            z = sExp(e,T,lnge,λ)
            𝚬 = 𝚬+(e*z)
            𝚬sq = 𝚬sq+((e^2)*z)
            𝚭 = 𝚭 + z
        end
        
        uT = (𝚬/𝚭) # Internal energy for the given temperature.
        cT = ((𝚬sq/𝚭)-(uT^2))/(T^2) 
        fT = -T*(λ+log(𝚭))
        entropyS = (uT-fT)/T
        
        us[i] = uT/n # Energy per spin.
        cs[i] = cT/n # Specific heat per spin.
        Fs[i] = fT/n
        Ss[i] = entropyS/n
    end

    writedlm(pathname*"/temperatures.csv",temperatures,',')
    writedlm(pathname*"/us.csv",us,',')
    writedlm(pathname*"/cs.csv",cs,',')
    writedlm(pathname*"/Fs.csv",Fs,',')
    writedlm(pathname*"/Ss.csv",Ss,',')

end






#thermo_WL(0.01,2.5,400,"simu6")


pathname = "./outputWL/simu6/"
temps = vec(readdlm(pathname*"/temperatures.csv"))
us = vec(readdlm(pathname*"/us.csv"))
cs = vec(readdlm(pathname*"/cs.csv"))
Fs = vec(readdlm(pathname*"/Fs.csv"))
Ss = vec(readdlm(pathname*"/Ss.csv"))


# plot(temps,us,xlabel="T",ylabel="",label="u(T)",lw=2,color="green",alpha=0.8)
# plot(temps,cs,label="c(T)",lw=2,color="blue",alpha=0.8)

# plot(temps,Fs,xlabel="T",ylabel="F(T)",label="",lw=2,color="red",alpha=0.8)

# plot(temps,Ss,xlabel="T",ylabel="S(T)",label="",lw=2,color="cyan",alpha=0.8)