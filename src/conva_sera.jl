module conva_sera

import HTTP
import CSV
import DataFrames
import Statistics
import Plots

using DataFrames
using Statistics
using Plots
using Distributions
using StatsBase


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
function get_pacf(data, max_lag; ifplot = false)
    
    k = max_lag
    data = get_acf(data,k,ifplot = false)
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


"""
        arma_solve(data,p,q)

# Arguments

- `data::Array`: The data want to be calculated.
   The data should be 1-Dimension Array.

- `p::Int`: the order of AR model.
   
- `q::Int`: the order of MA model.

# return 
        [a,b sigma]
- `a::Array`: coefficient of AR model
- `b::Array`: coefficient of MA model
- `sigma::Float`: sigma


# Examples
    

        output = arma_solve(data,a,b)
        a = output[1]
        b = output[2]
        sigma = output[3]

"""
function arma_solve(data,p,q)

    Gpq = zeros(p,p)
    pmax = round(sqrt(length(data)))

    gms = autocorvariance(data, trunc(Int,pmax+p))
    
    for i = 1:p
        for j = 1:p
            Gpq[i,j] = gms[abs(q+i-j)+1]
        end
    end

    gs = gms[(q+1+1):(q+p+1)]
    a = inv(Gpq)*gs
    aa = pushfirst!(a,-1)
    gys  = zeros(1,q+1)
        for k = 0:q
        sum1 = 0
        for i = 0:p
            for j = 0:p
                out = aa[i+1]*aa[j+1]*gms[abs(k+j-i)+1]
                sum1 = sum1 + out
            end
        end
        gys[k+1] = sum1
    end
    
    res = ma_solve(gys)
    b = res[1]
    sigma = sqrt(res[2])
    return [a[2:end],b[1:end],sigma]
end

"""
        autocorvariance(data, k)

# Arguments

- `data::Array`: The data want to be calculated.
   The data should be 1-Dimension Array.

- `k::Int`: max lag.

# return 

- `acf::Array`: the autocorvariance of the data.

# Examples

        data = [1,2,3,4,5,6,7,8,9,10]
        autocorvariance(data, 4)

"""
function autocorvariance(data,k)
    
    data = values(data)
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
        
        while i < n+1
            
            rk1 = (data[i]-avg) * (data[i-x]-avg)
            sum_rk1 = sum_rk1+rk1
            i = i+1
        end
        sum_rk1 = sum_rk1/n
       
        rk = sum_rk1
        rk = round(rk; digits =3 )  
        append!(acf, rk)  
    end
    return acf
end


"""
        ma_solve(gms)

# Arguments

- `gms::Array`: The autocorvariance of the data.
   The Array should be 1-Dimension Array.


# return 

[b sigma]

- `b::Array`: coefficient of MA model
- `sigma::Float`: sigma

# Examples
    
        output = ma_solve(gms)

        b = output[1]
        sigma = output[2]

"""
function ma_solve(gms)
    k = 100
    q = length(gms)-1
    if q == 1
        rho1 = gms[2]/gms[1]
        b = (1-sqrt(1-4*rho1^2))/(2*rho1)
        s2 = gms[1]/(1+b^2)
        return [b,s2]
    end
    
    A = zeros(q,q)
    
    for j = 2:q
        A[j-1,j] = 1
    end
    
    cc = zeros(q,1)
    cc[1] = 1
    gamma0 = gms[1]
    gammas = zeros(q+k)
    gammas[1:(q+1)] = gms
    gamq = gms[2:end]
    Gammak = zeros(k,k)
    
    for i = 1:k
        for j = 1:k
            Gammak[i,j] = gammas[abs(i-j)+1]
        end
    end
    
    Omk = zeros(q,k)
    for i = 1:q
        for j = 1:k
            Omk[i,j] = gammas[i+j-1+1]
        end
    end
    
    PI = Omk*(inv(Gammak)*transpose(Omk))
    s2 = gamma0 - (transpose(cc)*PI*cc)[1]
    b = 1/s2*(gamq-A*PI*cc)
    return [b,s2]
end

"""
        ma_gen(n,a;sigma=1)

# Arguments

- `n::Int`: the number of created data.
- `a::Array`: the coefficient of AR model.
- `sigma::Float`: sigma of model.

# return 

x

- `x::Array`: created data

# Examples
    
        output = ma_gen(10,a,sigma = 0.88)

"""
function ma_gen(n,a;sigma=1)
    q = length(a)
    n2 = n+q
    eps = rand(Normal(0,sigma),n2)
   
    w = pushfirst!(a,1.0)
    x2 = filter(eps,w,"convolution")
    x = x2[(q+1):n2]
    return x
end

function arma_gen(n,a,b;sigma = 1, n0=1000, x0 = zeros(length(n0),1))
    n2 = n0+n
    p = length(a)
    eps = ma_gen(n2,b,sigma = sigma)
    x2 = filter(eps,a,"recursive")
    x = x2[(n0+1):n2]
    return x
end

"""
        filter(x, w, method)

Applies linear filtering to a univariate time series or to each series separately of a multivariate time series.

# Arguments

- `x::Array`: one dimension data.
- `w::Array`: the coefficients.
- `method::String`: "convolution" or "recursive".

# return 

out

- `out::Array`: created data

# Examples
```julia
        x = [1,2,3,4,5,6]
        w = [1,2,3]

        output = filter(x,w,method = "convolution")
```
"""
function filter(x, w, method)
    lenx = length(x)
    lenw = length(w)
    if method == "convolution"
        out = zeros(lenx,1)
        for i =lenw:lenx
            for j = 1:lenw
                out[i] = out[i] + w[j]*x[i-j+1]
            end
        end
    end  
    
    if method == "recursive"
        lenx = lenx+1
        out = zeros(lenx,1)
        out[2] = x[1]
        if lenw>1
            for i = 1:lenw
                for j = 1:i
                    out[i+2] = w[j]*out[i+1-j+1]+out[i+2]
            
                end
                out[i+2] = out[i+2]+x[i+1]
            end
        end
        for i = lenw+1:lenx-2
            for j = 1:lenw
                out[i+2] = w[j]*out[i+1-j+1]+out[i+2]
            end
            out[i+2] =out[i+2]+x[i+1]
        end
        out = out[2:end]
    end
                        
    return out
end

"""
        arma_wold(n,a;b = zeros(1,1))

Wold coefficients for the ARMA model

# Arguments

- `n::Int`: number of data points.
- `a::Array`: the AR model coefficients.
- `b::Array`: the MA model coefficients.

# return 

psi

- `psi::Array`: Wold coefficients

# Examples

        a = [1.1, 0.2, 0.3]
        b = [0.2, -0.3]

        output = arma_wold(10,a,b = b)

"""
function arma_wold(n,a;b = zeros(1,1))
    p = length(a)
    q = length(b)
    arev = reverse(a)
    psi = zeros(n,1)
    psi[1] = 1
    for j = 1:n-1
        if j <= q
            bj = b[j]
        else
            bj = 0
        end
        psis = psi[max(1,j+1-p):j]
        np = length(psis)
        if np < p
            psis = append!(repeat([0.0],p-np),psis)
        end
        psi[j+1] = bj+sum(arev.*psis)
    end
    return psi
end

"""
        arma_gamma(n,a;b::Array = zeros(1,1),sigma = 1.0)

Calculate theoretical autocovariance function
of ARMA model using Wold expansion

# Arguments

- `n::Int`: number of data points/lags.
- `a::Array`: the AR model coefficients.
- `b::Array`: the MA model coefficients.
- `sigma::Float`: the sigma of the model

# return 

gam

- `gam::Array`: autocorvariance

# Examples

        a = [1.1, 0.2, 0.3]
        b = [0.2, -0.3]

        output = arma_gamma(10,a,b = b,sigma=0.95)

"""
function arma_gamma(n,a;b::Array = zeros(1,1),sigma = 1.0)
    nn = n + 100
    psi = arma_wold(nn,a;b=b)
    gam = zeros(n,1)
    for i = 0:n-1
        gam[i+1] = sum(psi[1:(nn-i)].*psi[(i+1):nn])
    end
    gam = (sigma^2)*gam
    return gam[1:end]
end


"""
        ts_ipred_coef(gam)

Using autocorvariance to calculate coefficients of best linear predictor

# Arguments

- `gam::Array`: the autocorvariance of the model.

# return 

[theta,nu]

- `theta::Array`: coefficients.
- `nu::Array`: square error of the predictor
# Examples

        gam = [1, 0.4, -0.4, 1.2, 2.3]
        output = ts_ipred_coef(gam)
        theta = output[1]
        nu = output[2]

"""
function ts_ipred_coef(gam)
    
    if ndims(gam) == 1
        n = length(gam)
        out = zeros(n,n)
        for i= 1:n
            for j = 1:n
                out[i,j] = gam[abs(i-j)+1]
            end
        end
        gam = out
    else
        n = size(gam)[1]
    end
    
    theta = zeros(n-1,n-1)
    nu = zeros(n,1)
    nu[1] = gam[1,1]
    theta[1,1] = gam[2,1]/gam[1,1]
    nu[2] = gam[2,2]-theta[1,1]^2*nu[1]
    for k = 2:n-1
        theta[k,k] = gam[k+1,0+1]/nu[1]
        for j = 1:k-1
            theta[k,k-j] = (gam[k+1,j+1]-sum(reverse(theta[j,1:j]).*reverse(theta[k,(k-j+1):k]).*nu[1:j]))/nu[j+1]
        end
        nu[k+1] = gam[k+1,k+1] - sum(reverse(theta[k,1:k].*theta[k,1:k]).*nu[1:k])
    end
    return [theta,nu]
end


"""
        ts_ipred(x,gams;demean = false, endx = length(x), conf = 0.95)

Using autocorvariance to do one-step-ahead prediction

# Arguments

- `x::Array`: data points.
- `gams::Array`: the autocorvariance of the model.
- `demean::Bool`: True if the mean is non-zero.
- `endx::Int`: the lenght of predictive values.
- `conf::Float`: confidence.

# return 

[x[1:end],pred[1:end],lb[1:end],ub[1:end],conf]

- `x::Array`: true data points.
- `pred::Array`: predictive values.
- `lb::Array`: the low boundary of predictive values
- `ub::Array`: the up boundary of predictive values
- `conf::Float`: the confidence.
```
        x = [1,2,3,4,5,6,7,8,9,10]
        gam = [1, 0.4, -0.4, 1.2, 2.3]
        output = ts_ipred(x,gams;demean = true, endx = length(x), conf = 0.95)

        x = output[1]
        pred = output[2]
```
"""
function ts_ipred(x,gams;demean = false, endx = length(x), conf = 0.95)
    if ndims(gams)==1
        stationary = true
    else
        stationary = false
    end
    res1 = ts_ipred_coef(gams)
    theta = res1[1]
    nu = res1[2]
    
    if demean == true
        xmean = mean(x)
        x = x.-xmean
    end
    
    if stationary == false
        endx = length(nu)
    end
    w = zeros(endx,1)
    pred = zeros(endx,1)
    lb = zeros(endx,1)
    ub = zeros(endx,1)
    lam = quantile.(Normal(0,1),1 - (1-conf)/2)
    
    pred[1] = 0
    w[1] = x[1]
    for i = 1:endx-1
        pred[i+1] =  sum(reverse(theta[i,1:i]).*w[1:i])
        w[i+1] = x[i+1] - pred[i+1]
    end
    sqrtnu = zeros(length(nu),1)
    for i = 1:length(nu)
        sqrtnu[i] = sqrt(nu[i])
    end
        
    lb = pred.-lam*sqrtnu
    ub = pred.+lam*sqrtnu
    
    if demean == true
        pred = pred.+xmean
        lb = lb.+xmean
        ub = ub.+xmean
    end
    
    return [x[1:end],pred[1:end],lb[1:end],ub[1:end],conf]  
end


"""
        arma_predict(x,a,b,sigma,n)

Using ARMA model to do prediction

# Arguments

- `x::Array`: data points.
- `a::Array`: coefficients of AR model
- `b::Array`: coefficients of MA model.
- `sigma::Float`: sigma of the ARMA model.
- `n::Int`: number of predictive values.

# return 

output

- `output::Array`: predictive data points.

# Example
```
        x = [1,2,3,4,5,6,7,8,9,10]
        a = [1.1, 0.4, -0.4]
        b = [1.2, 0.3]
        output = arma_predict(x, a, b, 1, 3)
```

"""
function arma_predict(x,a,b,sigma,n)
    
    lena = length(a)
    lenb = length(b)
    lenx = length(x)
    meanx = 0
    
    if mean(x)!= 0
        meanx = mean(x)
        x = x.-meanx
    end
    
    gams = arma_gamma(lenx,a,b = b,sigma = sigma)
    pred = ts_ipred(x,gams,demean= false,endx = lenx, conf = 0.95)
    pred = pred[2]
    residual = x.-pred
    append!(residual,zeros(n,1))
    append!(x,zeros(n,1))
    
    for i=1:n
        for q = 1:lena
            x[i+lenx] = x[i+lenx]+a[q]*x[i+lenx-q]
        end
        
        for p = 1:lenb
            x[i+lenx] = x[i+lenx]+b[p]*residual[i+lenx-p]
        end
    end
    output = x.+meanx
    return output[(lenx+1):(lenx+n)]
end


export get_acf,get_pacf, countryData, confData, deathsData, recData, arma_predict, ts_ipred_coef, ts_ipred, arma_gamma, arma_wold, filter, arma_gen, ma_gen
export ma_solve, arma_solve, autocorvariance

end # module

