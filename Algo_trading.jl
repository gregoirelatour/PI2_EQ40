#using QuandlAccess
#using Gadfly
using AlphaVantage
using DataFrames, Dates, Tables, TimeSeries, Statistics
using DataFramesMeta
using Plots
using Random
using PoissonRandom, Distributions
theme(:juno)

client = AlphaVantage.GLOBAL[]
client.key = "6065Z0SK19ZA64GE"
#6065Z0SK19ZA64GE api key alpha vantage

function raw_to_dataframe(rawData)
    df = DataFrame(rawData[1])
    dfNames = Symbol.(vcat(rawData[2]...))
    df = rename(df, dfNames)

    df.Date = Date.(df.timestamp)
    for x in (:open, :high, :low, :close, :adjusted_close, :dividend_amount)
        df[!, x] = Float64.(df[!, x])
    end 
    df.volume = Int64.(df.volume)
    return df
end

function intra_to_dataframe(rawData)
    df = DataFrame(rawData[1])
    dfNames = Symbol.(vcat(rawData[2]...))
    df = rename(df, dfNames)

    df.DateTime = DateTime.(df.timestamp, "yyyy-mm-dd HH:MM:SS")
    for x in (:open, :high, :low, :close)
        df[!, x] = Float64.(df[!, x])
    end 
    df.volume = Int64.(df.volume)
    return df
end

function returns_dataframe(data)
    y = Vector{Float64}(undef, length(data."open"))
    for i in range(1,length(data."open")-1)
        y[i]= (data."open"[i] - data."open"[i+1])/data."open"[i+1]
    end
    y[5000]= 0.0
    data."returns"= y
    
end

function cumulative_returns_dataframe(data)
    y = Vector{Float64}(undef, length(data."open"))
    y[5000]=0.0
    for i in range(1,length(data."open")-1)
        y[5000-i]= data."returns"[5000-i] + y[5001-i]
    end
    data."cumu_returns"= y
end

function mov_mean(data, min, max)
    res = Vector{Float64}(undef, length(data))
    for i in range(min + 1, length(data)-max)
        res[i]=mean(data[i-min: i+max])
    end
    
    return res
end


function C(m)
    return  2*log(m)^(1/2)-(log(pi)+log(log(m)))/(2*(2*log(m))^(1/2))
    
end
function S(m)
    return 1/2*(log(m)^(1/2))
    
end

function IV(t,returns,n)
    coef=(pi/2)*(n/(n-1))*(1/(n-1))
    res=0
    for i in range(n-2)
        res=res+abs(returns[t-i])*abs(returns[t-i-1])
    end
    return res*coef
    
end

function presence_jumps(returns,Am,m,n,eps)
    value = Vector{Float64}(undef, length(returns))
    Jump=Vector{Bool}(undef, length(returns))
    for i in range(Am)
        value[i]=(abs(returns[i]/IV(i,returns,n))-C(m))/S(m)
    
        Max=max(value[i])
        t=findall(value[i].=Max)
        psi=exp(exp(-t))
        if psi+eps<Max<=psi+eps
            Jump[i]=True
        else 
            Jump[i]=False
        end
    end
    return Jump
end

# function LeeMykland(data, delta, significance_level = 0.01)
#     tm = 252*24*60 #utile en 1 minute
#     k = ceil(sqrt(tm/delta))
#     r = append!(diff(log.(data)))
#     bpv = abs.(r) .* abs.(reverse(r))
#     bpv = 

    
# end

#Creation of the Option object


##CLUSTERING
function RandomPoisson(lambda, x, y)
    a = zeros(Int64, x,y)
    for i in range(1,x)
        a[i, :]=Array([pois_rand(lambda) for n in 1:y])
    end 
    return a 
end
function RandomNormal(mu,sigma, x, y)
    a = zeros(Float64, x,y)
    for i in range(1,x)
        a[i, :]=rand(Normal(mu, sigma), y)
    end 
    return a 
    
end
function GenererCourbesPoisson(nbTrajet,index,T,lp) #index coccorepond à l'increment, lp=paramettre, intensité lambda, T= temps, maturité
    #creation de matrices vides 
    #X=correspond au poisson et XC= au poisson compensé
    X= zeros(Float64, nbTrajet, index+1) #nombre de trajectoire en ligne et nombre l'increment en colonne
    XC= zeros(Float64, nbTrajet, index+1)
    temps = zeros(Float64, index+1)

    deltaT = T/index #mettre un float? correspond au pas

    P= RandomPoisson(lp*deltaT,nbTrajet,index)
    #moyenne 0 et variance 1. Variable ne sera plus biasé car on compense donc on aura une moyenne de 0
    for i in range(1,index)
        X[:,i+1]=X[:,i].+P[:,i]
        XC[:,i+1]=XC[:,i].-lp*deltaT .+ P[:,i] #on compense en soustrayant cf demos teams
        temps[i+1]=temps[i]+deltaT #on incremente de delta, on saute
    end
    trajectoires=Dict("temps"=>temps,"X"=>X,"XCompense"=>XC) #dictionnaire avec les differentes matrices 
    return trajectoires


end

function CalculMain()
    nbTrajet=28
    index=508
    
    T=38
    lp=1
    trajectoires=GenererCourbesPoisson(nbTrajet,index,T,lp)
    grilleTemps=trajectoires["temps"]
    X=trajectoires["X"]
    XC=trajectoires["XCompense"]

    # println(grilleTemps)
    # println(X)
    # println(XC)
    Plots.plot(grilleTemps,X')

    Plots.plot(grilleTemps,XC')

end

function GenererTrajetIto(nbTrajet,index,S0, T,lp,muJ,sigmaJ,r,sigma)
    #muJ= magnitude du saut
    #sigmaj= c'st lintertidue sur la taille du saut, la volatinité
    #sigma =volatité pour un mouvement brownien
    #S0=valeur initial du stock, prix depart
    #r= taux d'interet

    # creation matrice vide pour le processus de poisson et processus de poisson compensé
    X= zeros(Float64, nbTrajet, index+1)
    XC= zeros(Float64, nbTrajet, index+1)
    temps = zeros(Float64, index+1)

    deltaT = T / index #float ? #le pas de lincrement, du saut
    X[:,1] .= log(S0) # colonne 1, on met le logarithme du prix depart stock
    XC[:,1] .= S0

     # expérance E(e^J) pour J~N(muJ,sigmaJ^2)
    EeJ = exp(muJ + 0.5*sigmaJ*sigmaJ) #esperance d'une distribution log normal calculé analytiquement cf teams
    P= RandomPoisson(lp*deltaT,nbTrajet,index)

    Z = RandomNormal(0.0,1.0,nbTrajet,index)  # increment pour la distribution normal, donc ici mouvement brownin

    J = RandomNormal(muJ,sigmaJ,nbTrajet,index) #VA sur l'amplitude du saut, correspont à lespérance au dessus 
    for i in range(1,index)
        # on veut espérance nul et variance 1. donc on compense comme précédement 
        #on normalise
        if nbTrajet > 1
            Z[:,i] = (Z[:,i] .- mean(Z[:,i])) ./ std(Z[:,i]) # on centre (soustraction ) et on reduit (division)
        end
        X[:,i+1]  = X[:,i] .+ (r - lp*(EeJ-1) - 0.5*sigma*sigma)*deltaT +sigma*sqrt(deltaT) .* Z[:,i] .+ J[:,i] .* P[:,i] # formule analytique cf teams
        temps[i+1] = temps[i] + deltaT # on incremente, on saute
    end
    XC = exp.(X)
    trajectoires=Dict("temps"=>temps,"X"=>X,"XC"=>XC)
    return trajectoires
end

function CalculMain2()
    nbTrajet = 28
    index = 508
    T = 38
    lp = 1
    muJ = 0.3
    sigmaJ = 0.7
    sigma = 0.3

    S0 =100
    r=0.05
    trajectoires = GenererTrajetIto(nbTrajet,index,S0, T,lp,muJ,sigmaJ,r,sigma)
    grilleTemps = trajectoires["temps"]
    X = trajectoires["X"]
    XC = trajectoires["XC"]
    Plots.plot(grilleTemps,X')

    #Plots.plot(grilleTemps,XC')


end


mutable struct Option
    market
    type
    open #Days
    close
    returns

    function Option(market,type, open, close)
        this = new()
        this.type= type
        this.open = open
        this.close = close

        o = values(market[open])[1]
        c = values(market[close])[1]
        

        if type == "SHORT"
            this.returns = -(c/o-1)
        elseif type == "LONG"
            this.returns = c/o-1

        end

        return this
        
    end

end
# function pre_avrage_log_returns(data)
#     y = Vector{Float64}(undef, length(data."open"))
#     y[1]= 0.0
#     K= Int64(length(data."open")/2)
#     t=0.0
#     t-1 = 0.0
#     for i in range(2, length(data.open))
#         for j in range (2, i)
#             t = log(y[i-K+j])
    
# end
function testEURUSD()
    eurusdRaw = AlphaVantage.fx_daily("EUR", "USD",outputsize="full")
    EURUSD = DataFrame(eurusdRaw)
    EURUSD[!, :timestamp] = Dates.Date.(EURUSD[!, :timestamp])
    timearr = TimeArray(EURUSD, timestamp = :timestamp)

    data_open = timearr[:open]
    #x = data_open[Date(2002,12,16)] - data_open[Date(2002,12,20)]
    y = Vector{Float64}(undef, length(EURUSD[!, :open]))
    y = mov_mean(EURUSD[!, :open],0,10)

    short1 = Option(data_open,"SHORT", Date(2022,02,03),Date(2022,02,09))
    long1 = Option(data_open,"LONG", Date(2022,02,03),Date(2022,02,09))

    println(short1.returns," ", long1.returns)
    Plots.plot(y)
end


CalculMain2()
# returns_dataframe(EURUSD)
# cumulative_returns_dataframe(EURUSD)




#Plot
# Plots.plot(EURUSD[!, :timestamp], EURUSD[!, :open], label=["Open"], dpi=500,thickness_scaling=1, minorgrid = true)


###test avec l'api de QUANDL####
#quandl = Quandl("wTCUEtx5HUzJFV28QAQq")

#test=quandl(TimeSeries("ML/BBY"))

#spx = quandl(TimeSeries("BITFINEX/EUR/USD"))

#Gadfly.plot(x=spx.Date, y=spx.High, Geom.line)
