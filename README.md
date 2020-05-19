## How to install the package

```julia
using Pkg
Pkg.add(PackageSpec(url = "https://github.com/a365020121/Covid19Data.jl"))
using conva_sera
```

The functions and Dataset you can use now:
```julia
function:

countryData()
get_acf()
arma_predict()
ts_ipred_coef
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

How to use this funcitons:

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


```
