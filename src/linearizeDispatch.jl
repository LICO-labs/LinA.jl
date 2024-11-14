

#This file is used to dispatch Linearize() on an optional argument (e::ErrorType)

"""
    Linearize(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType; bounding = Best() ::BoundingType, ConcavityChanges = [Inf]::Array{Float64,1})

Makes an optimal piecewise Linear approximation of expr_fct from x1 to x2. The result will be an array of `LinearPiece`. Note that the array is directly callable as a function.
# Arguments
- `expr_fct` : function to linearize from R to R (either an expression or a native julia function)
- `x1` : from
- `x2` : to
- `e` : error type either Absolute() or Relative()

# Optional Arguments
- `bounding` : `Under()` for an underestimation, `Over()` for an overestimation, `Best()` for estimation that can go under or over the function. By default, it uses `Best()`
- `ConcavityChanges` : Concavity changes in the function. If not given, they will be computed automatically which, in rare cases, can lead to precision errors if the concavity is asymptotic to zero. 

!!! note
    It is also possible to specify which algorithm to use between `HeuristicLin()` and `ExactLin()` by simply adding it after the error type.
    By default LinA uses the heuristic.


!!! note
    If the function is given by a expression, the variable is assume to be `x`

# Example
```julia-repl
julia> pwl = Linearize(:(x^2),0,2,Absolute(0.1))
3-element Vector{LinA.LinearPiece}:
 0.894427190999916 x -0.1 from 0.0 to 0.894427190999916
 2.683281572999748 x -1.7000000000000006 from 0.894427190999916 to 1.7888543819998326
 4.736067977499794 x -5.372135954999589 from 1.7888543819998326 to 2.0

julia> pwl(1)
0.9832815729997475
```
"""
function Linearize(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType; bounding = Best() ::BoundingType,
ConcavityChanges = [Inf]::Array{Float64,1})

    (isfinite(x1) && isfinite(x2)) || throw(ArgumentError("Must be called on a finite interval"))

    return ScaledLinearize(expr_fct, x1, x2, e, HeuristicLin, bounding, ConcavityChanges)
end

function Linearize(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType,algorithm::HeuristicLin; bounding = Best() ::BoundingType,
    ConcavityChanges = [Inf]::Array{Float64,1})

    (isfinite(x1) && isfinite(x2)) || throw(ArgumentError("Must be called on a finite interval"))
    
    return ScaledLinearize(expr_fct, x1, x2, e, HeuristicLin, bounding, ConcavityChanges)
end

function Linearize(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType,algorithm::ExactLin; bounding = Best() ::BoundingType,
    ConcavityChanges = [Inf]::Array{Float64,1})

    (isfinite(x1) && isfinite(x2)) || throw(ArgumentError("Must be called on a finite interval"))

    return ScaledLinearize(expr_fct, x1, x2, e, ExactLin, bounding, ConcavityChanges)
end

# The function ScaleLinearize will 
# 1) invert a negative function to make it positive; 
# 2) scale it so it is defined in the interval [0, 1] and such that max(f(x): x in [0, 1]) = 1.
function ScaledLinearize(f::Ef, x1::Real, x2::Real, e::ErrorType, LinAlg::Union{Type{ExactLin}, Type{HeuristicLin}}, bounding::BoundingType, concavity_changes)::Vector{LinearPiece}
    if is_mostly_negative(f, x1, x2)
        invert = -1
        g = invert_function(f)
        new_bounding = bounding isa Under ? Over() : (bounding isa Over ? Under() : Best())
    else 
        invert = 1
        g = f
        new_bounding = bounding
    end
    s = get_scale(g, x1, x2)
    h = scale_function(g, s, x1, x2)
    newe = e isa Absolute ? Absolute(e.delta / s) : e

    # find roots of f and prevent running the main alg at them
    rts = IntervalRootFinding.roots(x -> h(x), interval(0, 1))
    breakpoints = [0.0, 1.0]
    for z in rts
        zmid = (z.region.bareinterval.lo + z.region.bareinterval.hi) / 2 
        if zmid > 1e-5
            push!(breakpoints, zmid - 1e-5)
        end
        if zmid < 1 - 1e-5
            push!(breakpoints, zmid + 1e-5)
        end
    end
    sort!(breakpoints)
    lps = LinearPiece[]
    for i in 1:length(breakpoints) - 1
        xp0, xpf = breakpoints[i], breakpoints[i + 1]
        if xpf - xp0 < EPS
        elseif xpf - xp0 < 1e-4
            lp = construct_constant_piece(h, xp0, xpf, new_bounding)
            push!(lps, lp)
        else
            newlps = LinAlg(h, xp0, xpf, newe; bounding = new_bounding, ConcavityChanges = deepcopy(concavity_changes))
            # newlps = remove_infeasibilities(newlps, h, new_bounding)
            res = [optimize(x -> h(x) - lp(x), lp.xMin, lp.xMax) for lp in newlps]
            res = [minimum(r) for r in res]
            # println("res = $(res)")
            append!(lps, newlps)
        end
    end
    newlps = [invert_linearpiece(scale_linearpiece(lp, s, x1, x2), invert) for lp in lps]
    res = [optimize(x -> f(x) - lp(x), lp.xMin, lp.xMax) for lp in newlps]
    res = [minimum(r) for r in res]
    println("hola")
    println("res = $(res)")
    return newlps
end