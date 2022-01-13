include("convexConcaveSplit.jl")
include("CorridorFromInfo.jl")
include("ConvexCorridor.jl")



function HeuristicLin(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType; bounding = Best() ::BoundingType,
                    ConcavityChanges = [Inf]::Array{Float64,1},m = HeuristicLin() ::AbstractAlgorithm )
    

    if x1 >= x2
        return Array{LinearPiece}(undef, 0)
    end

   
    if ConcavityChanges == [Inf]
        ConcavityChanges = decoupageConcavite(x1,x2,expr_fct)
    end
    

    ConcavityChanges = [x1;ConcavityChanges;x2]
    ConcavityChanges = sort(unique(ConcavityChanges)) # make sure that the bounds are there



    pwl = Array{LinearPiece}(undef, 0)
    temp =  fctMaker(expr_fct)
    f(x) = temp(x)
    

    for i in 1 : length(ConcavityChanges)-1


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

        corridor = corridorFromInfo(x1,x2,expr,e,bounds)
        tempCorridor = LinearizeConvex(corridor...)
        
        inverted ? push!(pwl, -tempCorridor...) : push!(pwl, tempCorridor...)

    end


    return pwl

end