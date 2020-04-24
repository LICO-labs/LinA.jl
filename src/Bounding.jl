
include("Linearize Dispatch.jl")
include("exact method.jl")

function LinearBounding(expr_fct::Expr,x1::Real,x2::Real, e::Absolute;
    ConcavityChanges = [Inf]::Array{Float64,1} )

under = Linearize(expr_fct,x1,x2, e; bounding = Under(),ConcavityChanges = ConcavityChanges )

return under, under + 2*e.delta

end
function LinearBounding(expr_fct::Expr,x1::Real,x2::Real, e::ErrorType;
    ConcavityChanges = [Inf]::Array{Float64,1} )
lower = Linearize(expr_fct,x1,x2, e; bounding = Under(),ConcavityChanges = ConcavityChanges )
upper = Linearize(expr_fct,x1,x2, e; bounding = Over(),ConcavityChanges = ConcavityChanges )
return lower,upper 


end

