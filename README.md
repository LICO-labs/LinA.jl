# LinA.jl
LinA is a piecewise linear approximation package which approximates univariate diffentiable functions with an optimal piecewise linear function in term of number of pieces given a error metric.
Both absolute and relative errors are implemented. It is also possible to add custom error type. LinA work expression as well as with native Julia functions and
Details of the algorithms can be found in [this paper.](https://hal.archives-ouvertes.fr/hal-03336003)


# Usage of the LinA package

All the useful documentation should be accessible through the `?` native Julia command to access documentation. (ex. do ?Linearize() to access the documentation of that function)
## Basic use
To Linearize $f(x) = x^2$ from -10 to 10 with a maximum absolute error of $2$ simply do
```julia
julia> using LinA

julia> pwl = Linearize(:(x^2),-10,10,Absolute(2))
5-element Vector{LinA.LinearPiece}:
 -16.000000000000004 x -62.00000000000003 from -10.0 to -6.0
 -8.0 x -14.0 from -6.0 to -2.0
 0.0 x + 2.0 from -2.0 to 2.0
 8.0 x -14.0 from 2.0 to 5.999999999999999
 15.999999999999996 x -61.99999999999998 from 5.999999999999999 to 10.0

```
**Note:** by default LinA uses the hybrid heuristic algorithm. To use the (slightly slower) exact algorithm simply add `ExactLin()` as an argument.

You can now call `pwl` as a julia function such as

```julia
julia> pwl(2)
4.470129472588532
```
But also as an array to get the individual linear segments such as
```julia 
julia> pwl[2]
-8.0 x -14.0 from -6.0 to -2.0

julia> pwl[2].xMax
-2.0
```
## Plotting
Pwl functions are compatible with Plots.jl. To plot a pwl function simply do
```julia
using Plots

plot(x->pwl(x),-10,10)

```
![alt text](https://i.imgur.com/7IHj3qp.png)

## Citing

If you have used our librarie and wish to cite our work (which we greatly encourage) use the referrence of [our paper.](https://hal.archives-ouvertes.fr/hal-03336003) Starring the _LinA.jl_ repository on GitHub is also appreciated.
