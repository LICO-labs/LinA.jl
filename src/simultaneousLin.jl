

function SimultaneousLin(expr_fct::Union{Expr, Function}, x1::Real, x2::Real, e::LinA.ErrorType,
    algorithm; ConcavityChanges = [Inf]::Array{Float64,1})
    
    simultaneousLin(expr_fct,x1,x2,e;algorithm=algorithm, ConcavityChanges = ConcavityChanges )
    
end




"""
    SimultaneousLin(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType; ConcavityChanges = [Inf]::Array{Float64,1} )

Makes a pair of optimal piecewise Linear underestimation and overestimation of expr_fct from x1 to x2 that share the same breakpoints.
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
function SimultaneousLin(expr_fct::Union{Expr, Function}, x1::Real, x2::Real, e::LinA.ErrorType;
        algorithm = LinA.HeuristicLin(), ConcavityChanges = [Inf]::Array{Float64,1})
    
    
    if typeof(e) == LinA.Absolute
        return LinearBounding(expr_fct, x1, x2, e,algorithm; ConcavityChanges)
    end
    
    
    
    #note : this is inneficient as it recomputes the whole linearization between after every new linear piece
    #TODO : rewrite this efficiently (which requires to do different things for exact and heuristic) 
    tol = 0.00001 
    l = Array{LinA.LinearPiece}(undef, 0)
    o = Array{LinA.LinearPiece}(undef, 0)
    
    while x2 - x1 > tol
    
        underPwl = Linearize(expr_fct, x1, x2, e, algorithm;
            bounding = LinA.Under(), ConcavityChanges)[1]
        overPwl =  Linearize(expr_fct, x1, x2, e, algorithm;
            bounding = LinA.Over(), ConcavityChanges)[1]

        newx1 = max(underPwl.xMax,overPwl.xMax)
        push!(l,LinA.LinearPiece(x1,newx1,underPwl.a,underPwl.b,underPwl.fct))
        push!(o,LinA.LinearPiece(x1,newx1,overPwl.a,overPwl.b,overPwl.fct))

        x1=newx1
    end

    
    l,o
end



