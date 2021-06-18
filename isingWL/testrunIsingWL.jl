using Base: Int64
using DelimitedFiles

# Here I test the WL algorithm. PLot the internal energy per spin.

using Plots
gr()








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
    enSpinIsingWL(us)
Given an array of simulated internal energies per spin `us`; plots the given data along with the analytical solution over the temperature
range T∈[0,5]. The transition should be located around T=2.3.
"""
function enSpinIsingWL(us)

    tempsplt=range(0,stop=5,length=201)
    usplt=Float64[enspin(temp,1) for temp in tempsplt]

    plt=scatter(range(0,stop=5,length=length(us)),us,ms=5,color="green",xlabel="T",ylabel="u(T)",label="WL",title="WL for Ising model",alpha=0.7,legend=:topleft)
    plot!(tempsplt,usplt,label="Analytic",color="gray",alpha=0.7,lw=2)

    return plt
end









"""
    heatSpinIsingWL(cs)
Given an array of simulated heat capacities per spin `us`; plots the given data along with the analytical solution over the temperature
range T∈[0,5]. The transition should be located around T=2.3.
"""
function heatSpinIsingWL(cs)

    tempsplt=range(0,stop=5,length=201)
    csplt=Float64[capspin(temp,1) for temp in tempsplt]

    plt=scatter(range(0,stop=5,length=length(cs)),abs.(cs),ms=5,color="green",xlabel="T",ylabel="c(T)",label="WL",title="WL for Ising model",alpha=0.7,legend=:topleft)
    plot!(tempsplt,csplt,label="Analytic",color="gray",alpha=0.7,lw=2)

    return plt
end










"""
    freeEnIsingWL(fs)
Given an array of simulated free energies per spin `fs`; plots the given data over the temperature
range T∈[0,5]. The transition should be located around T=2.3.
"""
function freeEnIsingWL(fs)

    plt=scatter(range(0,stop=5,length=length(fs)),fs,ms=5,color="green",xlabel="T",ylabel="f(T)",label="WL",title="WL for Ising model",alpha=0.7,legend=:topleft)

    return plt
end












"""
    entIsingWL(Ss)
Given an array of simulated entropies per spin `Ss`; plots the given data over the temperature
range T∈[0,5]. The transition should be located around T=2.3.
"""
function entIsingWL(Ss)

    plt=scatter(range(0,stop=5,length=length(Ss)),Ss,ms=5,color="green",xlabel="T",ylabel="s(T)",label="WL",title="WL for Ising model",alpha=0.7,legend=:topleft)

    return plt
end





#Load the data. The first column are the temperatures, the second are the energies, the third the specific heats, the fourth the free energies,
# and the third are the entropies.
dataWL=readdlm("/Users/pedroruiz/Desktop/Diego/PF/Data/WLN16")





