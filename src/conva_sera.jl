module conva_sera

import HTTP
import CSV
import DataFrames
import Statistics
import Plots

using DataFrames
using Statistics
using Plots



confirmedFile =  "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
deathsFile = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
recoveredFile = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"

makeDF(url) = CSV.read(IOBuffer(String(HTTP.request("GET",url).body)))

confData = makeDF(confirmedFile)
deathsData = makeDF(deathsFile)
recData = makeDF(recoveredFile)


function cleanNames(df)
    rename!(df,Symbol("Province/State")=>:ProvinceOrState)
    rename!(df,Symbol("Country/Region")=>:CountryOrRegion)
    
    numDays = size(df)[2]-4 
    #first four fields are :ProvinceOrState, :CountryOrRegion, :Lat, :Long
    #first date is 1/22/20
    for i in 1:numDays
        rename!(df,names(df)[4+i]=>"day$i")
    end
end

cleanNames(confData)
cleanNames(deathsData)
cleanNames(recData)

"""
        countryData(country, dfType; sumProvinces = true)

# Arguments

- `country::String`: the name of country you want to choice from the dataframe.
  The country name should in the data.

- `dfType::Symbol`: there are three kinds of dfType: :confirmed, :deaths, :recovered.
   :confirmed: data of diagnoses
   :deaths: data of deaths
   :recovered: data of Rehabilitation

- `sumProvinces::Bool`: if the function should return the summary data of the country.
   For instance, 'sumProvinces = false'. Default argument is false.

# Examples
    
countryData("Australia", confirmed）
countryData("Australia", deaths, sumProvinces = false)

"""
function countryData(country,dfType; sumProvinces = true)
    if dfType == :confirmed
        df = confData
    elseif dfType == :deaths
        df = deathsData
    elseif dfType == :recovered
        df = recData
    else
        println("error")
    end
    dataMatrix = convert(Array{Int},df[df.CountryOrRegion .== country,5:end])
    if size(dataMatrix)[1] == 0
        println("error data country not found")
    end
    if sumProvinces
        return sum(dataMatrix,dims=1)[1,:]
    else
        return dataMatrix[1,:]
    end
end

"""
        get_acf(data, max_lag; ifplot = false)

# Arguments

- `data::Array`: The data want to be calculated.
   The data should be 1-Dimension Array.

- `max_lag::Int`: the number of lags.
   
- `ifplot::Bool`: plot the array or nor.

# Examples
    
data = [1,2,3,4,5,6,7,8,9,10]
get_acf(data, 2)
get_acf(data, 2, ifplot = true)

"""
function get_acf(data,max_lag; ifplot= false)
    
   k = max_lag
   acf = []
   n = length(data)
   avg = mean(data) 
    
   for x in (0:k)   
       rk = 0
       rk1 = 0
       rk2 = 0
       sum_rk1 = 0
       sum_rk2 = 0
       i = 1+x
       j = 1
       
       while i < n+1
           rk1 = (data[i]-avg) * (data[i-x]-avg)
           sum_rk1 = sum_rk1+rk1
           i = i+1
       end
        
       sum_rk1 = sum_rk1/n  
        
       while j < n+1
           rk2 = (data[j]-avg)^2
           sum_rk2 = rk2+sum_rk2
           j=j+1       
       end
        
       sum_rk2 = sum_rk2/n    
       rk = sum_rk1/sum_rk2
       rk = round(rk; digits =3 )  
       append!(acf, rk)       
   end
   
   if ifplot == true
       return plot(acf[2:end], seriestype = :scatter)
       
   elseif ifplot == false      
       return acf[2:end]
   end    
end

"""
        get_pacf(data, max_lag; ifplot = false)

# Arguments

- `data::Array`: The data want to be calculated.
   The data should be 1-Dimension Array.

- `max_lag::Int`: the number of lags.
   
- `ifplot::Bool`: plot the array or nor.

# Examples
    
data = [1,2,3,4,5,6,7,8,9,10]
get_pacf(data, 2)
get_pacf(data, 2, ifplot = true)

"""
function get_pacf(data, max_lag, ifplot = false)
    
    k = max_lag
    data = get_acf(data,k,false)
    result = zeros((length(data),length(data)))
    result[1,1] = data[1]
    
    for i = 1:length(data)-1
        sum_1 = 0
        sum_2 = 0
        
        for j = 1:i
            sum_1 = sum_1 + result[i,j]*data[i+1-j]
            sum_2 = sum_2 + result[i,j]*data[j]
        end
            
        result[i+1,i+1] = (data[i+1]-sum_1)/(1-sum_2) 
        
        for j = 1:i
            result[i+1,j] = result[i,j] - result[i+1,i+1]*result[i,i-j+1]
        end
    end
    
    pacf = []
    
    for i= 1:k
        append!(pacf, round(result[i,i],digits=3))
    end
    
    if ifplot == true
        return plot(pacf, seriestype = :scatter)
        
    elseif ifplot == false  
        return pacf
    end      
end

export get_acf, countryData, confData, deathsData, recData

end # module

