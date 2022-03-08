
include("Linearize Dispatch.jl")
include("exact method.jl")

function LinearBounding(expr_fct::Ef,x1::Real,x2::Real, e::Absolute;
    ConcavityChanges = [Inf]::Array{Float64,1} )

under = Linearize(expr_fct,x1,x2, e; bounding = Under(),ConcavityChanges = ConcavityChanges )

return under, under + 2*e.delta

end

"""
    LinearBounding(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType; ConcavityChanges = [Inf]::Array{Float64,1} )

Makes an optimal piecewise Linear underestimation and overestimation of expr_fct from x1 to x2. For certain error types, it can saves a lot of overhead over calling Linearize two times.
# Arguments
- `expr_fct` : function to linearize (either an expresion or a native julia function)
- `x1` : from
- `x2` : to
- `e` : error type. Both Absolute() and Relative() are implemented.

# Optional Arguments
- `ConcavityChanges` : Concavity changes in the function. If not given, they will be computed automatically which, in rare cases, can lead to precision errors if the concavity is asymptotic to zero. 

"""
function LinearBounding(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType; ConcavityChanges = [Inf]::Array{Float64,1} )
    
lower = Linearize(expr_fct,x1,x2, e; bounding = Under(),ConcavityChanges = ConcavityChanges )
upper = Linearize(expr_fct,x1,x2, e; bounding = Over(),ConcavityChanges = ConcavityChanges )
return lower,upper 
    
end

