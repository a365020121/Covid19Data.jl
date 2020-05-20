## How to install the package

```julia
using Pkg
Pkg.add(PackageSpec(url = "https://github.com/a365020121/Covid19Data.jl"))
using conva_sera
```

The functions and Dataset you can use now:
```julia
functions:

countryData()
get_acf()
arma_predict()
ts_ipred_coef()
ts_ipred()
arma_gamma()
arma_wold()
filter()
arma_gen()
ma_gen()
ma_solve()
arma_solve()
autocorvariance()

Dataset:

confData # confirmed dataset
deathsData # deaths dataset
recData # recovered dataset
```

How to use these funcitons:

## you also can add ? before a funciton to find the document.
## like ?arma_gen()

```julia
# the data we use
data

20-element Array{Any,1}:
 -0.00351557541645414
  1.72462813688346
 -2.19592312526672
  1.29734468931968
 -0.661953072172693
 -1.33466852065149
  1.9297827432988
 -1.53239037841357
  0.831036732606825
  0.7279087375061
 -1.92215676302932
  0.428023229759807
  2.55058674699842
 -0.345451999285548
 -2.46858202605668
  0.56150616425268
  0.53850088453093
  0.187654544274226
  0.572737194799738
 -0.933389480023021
```

```julia
get_acf(data, 4)

4-element Array{Any,1}:
 -0.411
 -0.286
  0.293
 -0.219
```
```julia
get_pacf(data,4)

4-element Array{Any,1}:
 -0.411
 -0.547
 -0.196
 -0.489

```

```julia
model1 = arma_solve(data,4,2)

3-element Array{Any,1}:
  Any[-0.8264399894755466, -1.322984543103093, -0.6402583265659944, -0.5686458302955866]
  [0.36177593981733186, -0.4840362849481408]                                            
 1.0964763844324286

```

```julia
a = model1[1]
b = model1[2]
sigma = model1[3]
arima_gen(20,a,b,sigma = sigma)

20-element Array{Float64,1}:
 -3.3002306756247632  
 -1.6983822929728514  
  1.101689146219191   
  2.0096562736887185  
  1.741935726607993   
 -3.374351111809462   
 -2.1800677328353704  
  3.4003478767331625  
  4.2078675982045475  
 -2.4131319251074754  
 -4.454512572913576   
  0.4231263737292279  
  4.18755766868917    
 -0.026278600714865874
 -4.249290048174604   
  1.4090129073415762  
  1.0392789518013756  
  1.1015365744614676  
  1.3887528350788414  
 -5.646164621052959   
```

```julia
x = [1,2,3,4,5,6]
w = [1,2,3]
conva_sera.filter(x,w,"convolution")

6×1 Array{Float64,2}:
  0.0
  0.0
 10.0
 16.0
 22.0
 28.0
```

```julia
conva_sera.filter(x,w,"recursive")
6-element Array{Float64,1}:
   1.0
   3.0
   8.0
  21.0
  51.0
 123.0
```

```julia
gams = arma_gamma(20,a,b=b ,sigma = sigma)
# autocorvariance of the data(model)
20-element Array{Float64,1}:
  7.713000000115598  
 -1.5530000000439377 
 -5.424000000056405  
  2.4820000000557516 
  1.7329999999909427 
 -0.3600000000052417 
 -0.49999999998731953
 -1.6314522005069965 
  1.2548193846216729 
  1.646194787246017  
 -1.6917140819773706 
 -0.6554801585538355 
  1.0122888666196947 
  0.1776263433351677 
 -0.10437724882823471
 -0.42412569099732844
 -0.20075666993997893
  0.6928469968927747 
  0.02390519893038773
 -0.566668644727381  
```

```julia
pred = arma_predict2(data,a,b,sigma,10)

10-element Array{Float64,1}:
 -1.0649558370926162 
  1.4419767303674473 
 -0.4743907577349372 
  1.1434238205068454 
 -0.6579499984832787 
 -1.5081528293967468 
  1.6315936349742772 
  0.39496787656242316
 -1.1681780335612468 
  0.23292049751688373

```
