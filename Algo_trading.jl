#using QuandlAccess
#using Gadfly
using AlphaVantage
using DataFrames, Dates
using DataFramesMeta
using Plots
theme(:juno)



client = AlphaVantage.GLOBAL[]
client.key = "6065Z0SK19ZA64GE"
#6065Z0SK19ZA64GE api key alpha

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

# tslaRaw = AlphaVantage.time_series_daily("", outputsize="full")

# Plots.plot(tslaRaw)

# gr(size=(800,470))
# spy = time_series_daily("SPY", outputsize="full");
# # Convert to a DataFrame
# data = DataFrame(spy);
# # Convert timestamp column to Date type
# data[!, :timestamp] = Dates.Date.(data[!, :timestamp]);
# data[!, :open] = Float64.(data[!, :open])
# # Plot the timeseries
# Plots.plot(data[!, :timestamp], data[!, :open], label=["Open"])
function returns_dataframe(data)
    y = Vector{Float64}(undef, length(data."open"))
    y[1]= 0.0
    for i in range(2, length(data."open"))
        y[i]= (data."open"[i] - data."open"[i-1])/data."open"[i-1]
    end
    data."returns"= y
    
end

function cumulative_returns_dataframe(data)
    y = Vector{Float64}(undef, length(data."open"))
    y[1]= 0.0
    for i in range(2, length(data."open"))
        y[i]= data."returns"[i] + y[i-1]
    end
    data."cumu_returns"= y
    
end
eurusdRaw = AlphaVantage.fx_daily("EUR", "USD",outputsize="full")
EURUSD = DataFrame(eurusdRaw)
EURUSD[!, :timestamp] = Dates.Date.(EURUSD[!, :timestamp])

returns_dataframe(EURUSD)
cumulative_returns_dataframe(EURUSD)
Plots.plot(EURUSD[!, :timestamp], EURUSD[!, :open], label=["Open"], dpi=1000,thickness_scaling=1, minorgrid = true)


###QUANDL####
#quandl = Quandl("wTCUEtx5HUzJFV28QAQq")

#test=quandl(TimeSeries("ML/BBY"))

#spx = quandl(TimeSeries("BITFINEX/EUR/USD"))

#Gadfly.plot(x=spx.Date, y=spx.High, Geom.line)