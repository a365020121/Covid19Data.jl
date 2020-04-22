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

Dataset:

confData # confirmed dataset
deathsData # deaths dataset
recData # recovered dataset
```
