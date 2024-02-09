



#to dispatch
function LinearBounding(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType;
    ConcavityChanges = [Inf]::Array{Float64,1} )
    
    LinearBounding(expr_fct,x1,x2, e , HeuristicLin());ConcavityChanges = ConcavityChanges  )
    
end


function LinearBounding(expr_fct::Ef,x1::Real,x2::Real, e::Absolute,algorithm;
    ConcavityChanges = [Inf]::Array{Float64,1} )

under = Linearize(expr_fct,x1,x2, e; bounding = Under(),algorithm,ConcavityChanges = ConcavityChanges )

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

!!! note
    It is also possible to specify which algorithm to use between `HeuristicLin()` and `ExactLin()` by simply adding it after the error type.
    By default LinA uses the heuristic.

"""
function LinearBounding(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType,algorithm; ConcavityChanges = [Inf]::Array{Float64,1} )
    
lower = Linearize(expr_fct,x1,x2, e; bounding = Under(),algorithm,ConcavityChanges = ConcavityChanges )
upper = Linearize(expr_fct,x1,x2, e; bounding = Over(),algorithm,ConcavityChanges = ConcavityChanges )
return lower,upper 
    
end

