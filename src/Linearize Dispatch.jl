include("Heuristic.jl")
include("exact method.jl")


function Linearize(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType; bounding = Best() ::BoundingType,
ConcavityChanges = [Inf]::Array{Float64,1})

    return HeuristicLin(expr_fct,x1,x2, e; bounding = bounding, ConcavityChanges = ConcavityChanges)
end

function Linearize(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType,algorithm::HeuristicLin; bounding = Best() ::BoundingType,
    ConcavityChanges = [Inf]::Array{Float64,1})

        return HeuristicLin(expr_fct,x1,x2, e; bounding = bounding, ConcavityChanges = ConcavityChanges)
end

function Linearize(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType,algorithm::ExactLin; bounding = Best() ::BoundingType,
    ConcavityChanges = [Inf]::Array{Float64,1})

        return exactLin(expr_fct,x1,x2, e; bounding = bounding, ConcavityChanges = ConcavityChanges)
end
