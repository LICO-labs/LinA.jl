
"""
isContinuous(pwl, ε = 1e-5)

Determine whether a pwl function is continuous up to a numerical precision of ε.
# Arguments
- `plw` : pwl function
- `ε` : numerical precision used to detect if the end endpoints of two segments are equal (to detect discontinuities)

"""
function isContinuous(pwl, ε = 1e-5) 

    
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
function breakpoints(pwl, ε = 1e-5) 

    
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
    return max(abs(f(x1)), abs(f(x2)))
end

function scale_function(f::Ef, s::Real)::Ef
    if f isa Expr
        return :( eval(f) * s )
    elseif f isa Function
        return x -> s * f(x)
    end
end

function scale_linearpiece(lp::LinearPiece, s::Real)::LinearPiece
    return LinearPiece(lp.xMin, lp.xMax, s * lp.a, s * lp.b, x -> s * (lp.a * x + lp.b))
end