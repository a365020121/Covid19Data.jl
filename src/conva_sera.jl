module conva_sera

import HTTP
import CSV
import DataFrames


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


end # module

