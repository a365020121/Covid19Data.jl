module conva_sera

import HTTP
import CSV
import DataFrames



confirmedFile =  "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
deathsFile = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
recoveredFile = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"


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

function countryData(country,dfType; sumProvinces = true)
    if dfType == :confirmed
        df = confirmedDF
    elseif dfType == :deaths
        df = deathsDF
    elseif dfType == :recovered
        df = recoveredDF
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

export cleanNames, countryData, confirmedFile, deathsFile, recoveredFile

end # module

