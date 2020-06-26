include("convexConcaveSplit.jl")
include("CorridorFromInfo.jl")
include("ConvexCorridor.jl")

#using TimerOutputs
#const tot = TimerOutput()

function HeuristicLin(expr_fct::Expr,x1::Real,x2::Real, e::ErrorType; bounding = Best() ::BoundingType, 
                    ConcavityChanges = [Inf]::Array{Float64,1},m = HeuristicLin() ::AbstractAlgorithm )
    #reset_timer!(tot)
    
    if x1 >= x2
        return Array{LinearPiece}(undef, 0)
    end
    
   # @timeit tot "concave split " begin
        if ConcavityChanges == [Inf]
            ConcavityChanges = decoupageConvavite(x1,x2,expr_fct)
        end
    #end
    
    ConcavityChanges = [x1;ConcavityChanges;x2] 
    ConcavityChanges = sort(unique(ConcavityChanges)) # make sure that the bounds are there
    
    
    
    pwl = Array{LinearPiece}(undef, 0)
    temp =  mk_function(:((x -> $expr_fct)))
    f(x) = temp(x) 
    #@eval f(x) = $expr_fct

    for i in 1 : length(ConcavityChanges)-1
        
     #   @timeit tot "check concave " begin
            #these variables are specific to each sub-corridors
            x1 = ConcavityChanges[i]
            x2 = ConcavityChanges[i+1]
            expr = expr_fct
            bounds = bounding
        
            # if the function is concave, multiply by *-1 and switch under by over
            inverted = false
            if f''((x1+x2)/2) <0
                expr = -expr
                inverted = true
                bounds = -bounds
            end

     #   end
        
     #   @timeit tot "trouve pour convexe " begin
            corridor = corridorFromInfo(x1,x2,expr,e,bounds)
            tempCorridor = LinearizeConvex(corridor...)
     #   end
        
     #   @timeit tot "inversion" begin
           inverted ? push!(pwl, -tempCorridor...) : push!(pwl, tempCorridor...)
     #   end

    end

   # return pwl,tot
    return pwl
    
end



#function LinearBounding(expr_fct::Expr,x1::Real,x2::Real, e::Absolute;
#                    ConcavityChanges = [Inf]::Array{Float64,1} )
#    
#    under = Linearize(expr_fct,x1,x2, e; bounding = Under(),ConcavityChanges = ConcavityChanges )
#        
#    return under, under + 2*e.delta
#    
#end
#function LinearBounding(expr_fct::Expr,x1::Real,x2::Real, e::ErrorType;
#                    ConcavityChanges = [Inf]::Array{Float64,1} )
#    lower = Linearize(expr_fct,x1,x2, e; bounding = Under(),ConcavityChanges = ConcavityChanges )
#    upper = Linearize(expr_fct,x1,x2, e; bounding = Over(),ConcavityChanges = ConcavityChanges )
#    return lower,upper 
#    
#    
#end





