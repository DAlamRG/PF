using StatsBase: minimum, maximum
using Base: Int64
using DelimitedFiles
using Plots
gr()

# Here I test the WL algorithm. Plot the internal energy per spin.
include("./isingWL.jl")









"""
    secondElliptic(κ,ϕ)
"""
function secondElliptic(κ,ϕ) 
    return sqrt(1-(κ*sin(ϕ))^2)
end







"""
    firstElliptic(κ,ϕ)
"""
function firstElliptic(κ,ϕ)
    return 1/sqrt(1-(κ*sin(ϕ))^2)
end







"""
    centralInt(f,a,b,n)
Given a function `f`, an interval `[a,b]`, and a number of subintervals `n`; returns the numeric integral of `f` over the specified domain. 
"""
function centralInt(f,a,b,n)
    xs=range(a, stop=b, length=n+1)
    integral=0
    for i in 1:n
        el=(xs[i]+xs[i+1])/2
        integral=integral+f(el)
    end
    h=(b-a)/n
    integral=integral*h
    return integral
end








"""
    enspin(T,J)
Given a temperature `T` and a bonding constant `J`; returns the analytical solution for the internal energy per spin of the 2D ising model.
"""
function enspin(T,J)
    β=1/T
    u=2*β*J # Define useful parameterer.
    κ=2*((sinh(u))/((cosh(u))^2)) # Define useful parameterer.
    K1=centralInt(x->firstElliptic(κ,x),0,pi/2,100) # Define first elliptic integral.
    E1=centralInt(x->secondElliptic(κ,x),0,pi/2,100) # Define second elliptic integral.
    
    t1=-(2*J)*tanh(u)
    t2=-J*(((sinh(u))^2-1)/(sinh(u)*cosh(u)))
    t3=(2/pi)*K1-1
    total=t1+t2*t3
    return total
end









"""
    capspin(T,J)
Given a temperature `T` and a bonding constant `J`; returns the analytical solution for the heat capacity per spin of the 2D ising model.    
"""
function capspin(T,J)
    β=1/T
    u=2*β*J # DEfino éste param útil
    κ=2*((sinh(u))/((cosh(u))^2)) # Defino éste param útil
    K1=centralInt(x->firstElliptic(κ,x),0,pi/2,100) # Defino la primer integral elíptica
    E1=centralInt(x->secondElliptic(κ,x),0,pi/2,100) # Defino la segunda int el
    
    t1=(4*κ/(pi))
    t2=(β*J*coth(u))^2
    t3=1-(tanh(u))^2
    t4=pi/2+(2*(tanh(u))^2-1)*K1
    total=t1*t2*(K1-E1-t3*t4)
    return total
end













"""
    enSpinIsingWL(us,temps)
Given an array of simulated internal energies per spin `us` and the corresponding temperatures `temps`; plots the given data 
along with the analytical solution over the temperaturerange T∈[0,5]. The transition should be located around T=2.3.
"""
function enSpinIsingWL(us,temps)

    tempsplt=range(0,stop=5,length=201)
    usplt=Float64[enspin(temp,1) for temp in tempsplt]

    plt=scatter(temps,us,ms=5,color="green",xlabel="T",ylabel="u(T)",label="WL",title="WL for Ising model",alpha=0.7,legend=:topleft)
    plot!(tempsplt,usplt,label="Analytic",color="gray",alpha=0.7,lw=2)

    return plt
end









"""
    heatSpinIsingWL(cs,temps)
Given an array of simulated heat capacities per spin `us` and the corresponding temperatures `temps`; plots the given data 
along with the analytical solution over the temperature range T∈[0,5]. The transition should be located around T=2.3.
"""
function heatSpinIsingWL(cs,temps)

    tempsplt=range(0,stop=5,length=201)
    csplt=Float64[capspin(temp,1) for temp in tempsplt]

    plt=scatter(temps,abs.(cs),ms=5,color="green",xlabel="T",ylabel="c(T)",label="WL",title="WL for Ising model",alpha=0.7,legend=:topleft)
    plot!(tempsplt,csplt,label="Analytic",color="gray",alpha=0.7,lw=2)

    return plt
end










"""
    freeEnIsingWL(fs,temps)
Given an array of simulated free energies per spin `fs` and the corresponding temperatures `temps`; plots the given data over the temperature
range T∈[0,5]. The transition should be located around T=2.3.
"""
function freeEnIsingWL(fs,temps)

    plt=scatter(temps,fs,ms=5,color="green",xlabel="T",ylabel="f(T)",label="WL",title="WL for Ising model",alpha=0.7,legend=:bottomleft)
    plot!(temps,fs,lw=3,color="lightseagreen",label="",alpha=0.4)
    return plt
end












"""
    entIsingWL(Ss,temps)
Given an array of simulated entropies per spin `Ss` and the corresponding temperatures `temps`; plots the given data over the temperature
range T∈[0,5]. The transition should be located around T=2.3.
"""
function entIsingWL(Ss,temps)

    plt=scatter(temps,Ss,ms=5,color="green",xlabel="T",ylabel="s(T)",label="WL",title="WL for Ising model",alpha=0.7,legend=:topleft)
    plot!(temps,Ss,lw=3,color="lightseagreen",label="",alpha=0.4)

    return plt
end
















"""
    ising2D_thermo(lngE,ti,tf,l,N)
 Given a matrix containing the enegy densities `lngE`, an initial (final) temperature `ti` (`tf`), a length for the temperature array `l`,
 and a lattice size `N`; returns a matrix containing the thermodynamic variables´ values at each desired temperature.    
"""
function ising2D_thermo(lngE,ti,tf,l,N)

    temps=range(ti,stop=tf,length=l) 
    us=ones(Float64,l) # This array will contain the internal energy per spin for each of the temperatures in `temps`.
    cs=ones(Float64,l) # 
    Fs=ones(Float64,l)
    Ss=ones(Float64,l)

    lngEaux=Dict{Int64,Float64}(lngE[k,1] => lngE[k,2] for k in 1:size(lngE)[1])
    kk=keys(lngEaux)
    mink=minimum(kk)
    maxk=maximum(kk)
    mrest=min(lngEaux[mink],lngEaux[maxk])
    
    lngE=Dict{Int64,Float64}(lngE[k,1] => (lngE[k,2]-mrest) for k in 1:size(lngE)[1])

    for i in 1:length(temps)
        T=temps[i]
        𝚬=0 # Average energy.
        𝚬sq=0 # Average of the squared energy.
        𝚭=0 # Partition function, normalizes the average energy.
        λ=determineλ(lngE,T)
        for element in lngE
            e,lnge=element # Extract the energy and natural logarithm of the energy density.
            z=sExp(e,T,lnge,λ)
            𝚬=𝚬+(e*z)
            𝚬sq=𝚬sq+((e^2)*z)
            𝚭=𝚭+z
        end
        
        uT=(𝚬/𝚭) # Internal energy for the given temperature.
        cT=((𝚬sq/𝚭)-(uT^2))/(T^2) 
        fT=-T*(λ+log(𝚭))
        entropyS=(uT-fT)/T
        
        us[i]=uT/(N^2) # Energy per spin.
        cs[i]=cT/(N^2)
        # Specific heat per spin.
        Fs[i]=fT/(N^2)
        Ss[i]=entropyS/(N^2)
    end

    infoM=ones(Float64,(l,5))
    infoM[:,1]=temps
    infoM[:,2]=us
    infoM[:,3]=cs
    infoM[:,4]=Fs
    infoM[:,5]=Ss

    return infoM
end





















#Load the data. The first column are the temperatures, the second are the energies, the third the specific heats, the fourth the free energies,
# and the fifth are the entropies.
lngE=readdlm("./Data/WLlngE16")

dataWL=ising2D_thermo(lngE,0.1,5,52,16)

# include("testrunIsingWL.jl")
# enSpinIsingWL(dataWL[:,2],dataWL[:,1])
# heatSpinIsingWL(dataWL[:,3],dataWL[:,1])
# freeEnIsingWL(dataWL[:,4],dataWL[:,1])
# entIsingWL(dataWL[:,5],dataWL[:,1])




