
include("exactPiece.jl")

function exactLin(expr_fct::Expr,x1::Real,x2::Real, e::ErrorType; bounding = Best() ::BoundingType, 
                    ConcavityChanges = [Inf]::Array{Float64,1} )
    if x1 >= x2
        return Float64[]
    end
   
    
    if ConcavityChanges[1] == Inf
        ConcavityChanges = decoupageConvavite(x1,x2,expr_fct)
    end
    ConcavityChanges = [x1;ConcavityChanges;x2] 
    ConcavityChanges = sort(unique(ConcavityChanges)) # make sure that the bounds are there
    pwl = Array{LinearPiece}(undef, 0)
    cor = corridorFromInfo(x1, x2,expr_fct, e, bounding)[3:end-1]
    
    x2Temp = -1
    i=1
    while x1 < x2
        i = searchsortedfirst(ConcavityChanges,x1)
        ConcavityChanges[i] == x1 ? x2Temp = ConcavityChanges[i+1] : x2Temp = ConcavityChanges[i]
        #println(x1,"  ",x2Temp)
        pwl = [pwl;Linearize(expr_fct, x1, x2Temp, e,bounding = bounding,ConcavityChanges=[] )] 
        pwl[end].xMax == x2 && break
        #println(pwl)
        #println("hello") 
        pwl[end] = exactPiece(pwl[end].xMin,x2,cor...)
        pwl[end].xMax == x2 && break
        #println(pwl,"\n\n\n")
        x1 = pwl[end].xMax
    end
    return pwl
end