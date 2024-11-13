

"""
    LinearizeConvex(x1,x2,lower::Function,upper::Function,du::Function)

Makes an optimal piecewise Linear approximation from x1 to x2 of a convex corridor
# Arguments
- lower : function at the bottom of the corridor
- upper : function on top of the corridor
- du : derivative of the upper function

"""
function LinearizeConvex(x1,x2,e::ErrorType,lower::Function,upper::Function,du::Function)

    
    slope = 0.0
    b = 0.0
    pwl = Array{LinearPiece}(undef, 0)
    
    #Δ represent the distance between the linear function and the bottom of the corridor
    Δ(x, (slope, b)) = slope * x + b - lower(x)

    #distance between the linear function and the bottom of the corridor at the start of the coridor if tangeant in "a"
    f(a, x1) = upper(a) + du(a) * (x1 - a) - lower(x1)

    while x2 - x1 > EPS
        if f(x2, x1) >= 0.0
            slope = (upper(x2) - lower(x1)) / (x2 - x1)
            b = lower(x1) - slope * x1
            push!(pwl,LinearPiece(x1,x2,slope,b,ParaToFncLin(slope,b)))
            break
        end
    
        a = find_zero(x -> f(x, x1), x1, x2)
    
        # @assert(!isnan(a), "a is still NAN!, f(x1) = $(f(x1)), f(x2) = $(f(x2))")

        slope = du(a)
        b = lower(x1) - slope * x1

        #find the second end of the segment
        if Δ(x2, (slope, b)) <= 0.0
            nextX1 = find_zero(x -> Δ(x, (slope, b)), a, x2)
        else
            nextX1 = x2
        end

        push!(pwl,LinearPiece(x1,nextX1,slope,b,ParaToFncLin(slope,b)))
        x1 = nextX1
    end
    return pwl

end

function perform_binary_sarch_lowf(f, x1, x2)
    eps = 1e-3
    xp = x1
    for _ in 1:5
        t = 1 + eps
        while f(xp) * f(x2) > 0 && abs(f(xp)) < EPS
            xp = min(xp + t * EPS, x2)
            t *= (1 + eps)
        end
        eps *= 0.5
        t = 1 + eps
        while f(xp) * f(x2) < 0 || abs(f(xp)) > EPS
            xp = max(xp - t * EPS, x1)
            t *= (1 + eps)
        end
        eps *= 0.5
    end
    t = 1 + eps
    while f(xp) * f(x2) > 0 && abs(f(xp)) < EPS
        xp = min(xp + t * EPS, x2)
        t *= (1 + eps)
    end
    return xp
end
