#using QuandlAccess
#using Gadfly
using AlphaVantage
using DataFrames, Dates, Tables, TimeSeries, Statistics
using DataFramesMeta
using Plots
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

# function LeeMykland(data, delta, significance_level = 0.01)
#     tm = 252*24*60 #utile en 1 minute
#     k = ceil(sqrt(tm/delta))
#     r = append!(diff(log.(data)))
#     bpv = abs.(r) .* abs.(reverse(r))
#     bpv = 

    
# end

#Creation of the Option object
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
# returns_dataframe(EURUSD)
# cumulative_returns_dataframe(EURUSD)




#Plot
# Plots.plot(EURUSD[!, :timestamp], EURUSD[!, :open], label=["Open"], dpi=500,thickness_scaling=1, minorgrid = true)


###test avec l'api de QUANDL####
#quandl = Quandl("wTCUEtx5HUzJFV28QAQq")

#test=quandl(TimeSeries("ML/BBY"))

#spx = quandl(TimeSeries("BITFINEX/EUR/USD"))

#Gadfly.plot(x=spx.Date, y=spx.High, Geom.line)
