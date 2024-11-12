
"""
isContinuous(pwl, ε = 1e-5)

Determine whether a pwl function is continuous up to a numerical precision of ε.
# Arguments
- `plw` : pwl function
- `ε` : numerical precision used to detect if the end endpoints of two segments are equal (to detect discontinuities)

"""
function isContinuous(pwl, ε = EPS) 

    
    for i in 1:length(pwl)-1
        
        temp = pwl[i].xMax
        
        if pwl[i](temp) - pwl[i+1](temp) > ε
            return false
        end
        
    end
    
    true
    
end





"""
breakpoints(pwl, ε = 1e-5)

Outputs the breakpoints needed for several solvers (CPLEX, Gurobi,...) to natively model PWL functions. 
# Arguments
- `plw` : pwl function
- `ε` : numerical precision used to detect if the end endpoints of two segments are equal (to detect discontinuities)

"""
function breakpoints(pwl, ε = EPS) 

    
    bpx = [pwl[1].xMin]

    
    for i in 1:length(pwl)-1
        
        temp = pwl[i].xMax
        push!(bpx,temp)
        
        if pwl[i](temp) - pwl[i+1](temp) > ε
            push!(bpx , temp)
        end
        
    end
    
    push!(bpx,pwl[end].xMax)
    
    return bpx
    
    
end

function get_scale(f::Ef, x1::Real, x2::Real)::Float64
    resmax = optimize(x -> -f(x), x1, x2)
    resmin = optimize(x -> f(x), x1, x2)
    vmax = minimum(resmax)
    vmin = minimum(resmin)
    return max(abs(vmax), abs(vmin))
end

function scale_function(f::Ef, s::Real, x1::Real, x2::Real)::Ef
    if f isa Expr
        return :( eval(f) * s )
    elseif f isa Function
        return y -> f(x1 + y * (x2 - x1)) / s
    end
end

function scale_linearpiece(lp::LinearPiece, s::Real, x1::Real, x2::Real)::LinearPiece
    ymin = x1 + lp.xMin * (x2 - x1) 
    ymax = x1 + lp.xMax * (x2 - x1)
    ap = s * lp.a / (x2 - x1)
    bp = s * (lp.b - lp.a * x1 / (x2 - x1))
    return LinearPiece(ymin, ymax, ap, bp, x -> ap * x + bp)
end

function find_zeros(f::Ef, x1::Real, x2::Real)::Vector{Real}
    zs = IntervalRootFinding.roots(f, interval(x1, x2))
    if isempty(zs) return [] end
    zmin = argmin(z -> min(abs(f(z.region.bareinterval.lo)), abs(f(z.region.bareinterval.hi))), zs)
    return abs(f(zmin.region.bareinterval.lo)) < abs(f(zmin.region.bareinterval.hi)) ? [zmin.region.bareinterval.lo] : [zmin.region.bareinterval.hi]
end

function find_zero(f::Ef, x1::Real, x2::Real)::Real
    zeros = LinA.find_zeros(f, x1, x2)
    return isempty(zeros) ? NaN : zeros[begin]
end

function remove_infeasibilities(lps, g, bounding)
    newlps = []
    for lp in lps
        if bounding isa Under
            res = optimize(x -> g(x) - lp(x), lp.xMin, lp.xMax)
            d = minimum(res)
        elseif bounding isa Over
            res = optimize(x -> lp(x) - g(x), lp.xMin, lp.xMax)
            d = - minimum(res)
        else d = 0.0
        end
        if abs(d) > EPS
            push!(newlps, LinearPiece(lp.xMin, lp.xMax, lp.a, lp.b + d, x -> lp.a * x + lp.b + d))
        else
            push!(newlps, lp)
        end
    end
    return newlps
end