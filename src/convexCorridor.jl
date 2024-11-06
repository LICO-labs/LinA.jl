

"""
    LinearizeConvex(x1,x2,lower::Function,upper::Function,du::Function)

Makes an optimal piecewise Linear approximation from x1 to x2 of a convex corridor
# Arguments
- lower : function at the bottom of the corridor
- upper : function on top of the corridor
- du : derivative of the upper function

"""
function LinearizeConvex(x1,x2,lower::Function,upper::Function,du::Function)

    
    slope = 0.0
    b = 0.0
    pwl = Array{LinearPiece}(undef, 0)
    

    lin(x) = slope * x + b
    #Δ represent the distance between the linear function and the bottom of the corridor
    Δ(x) = lin(x) - lower(x)
    #distance between the linear function and the bottom of the corridor at the start of the coridor if tangeant in "a"
    f(a) = upper(a) + du(a)*(x1-a) - lower(x1)
    
    # println(stderr, "loop")
    while x2 - x1 > EPS
        # println(stderr, "asdasd")
        # println("entering loop with x1 = $x1, x2 = $x2, f(x1) = $(f(x1)), f(x2) = $(f(x2))")
        # if f(x2) < 0.0 && f(x1) * f(x2) > 0.0 break end
        #find the tangence point if possible
        if f(x2) >= 0.0
            slope = (upper(x2) - lower(x1)) / (x2 - x1)
            b = lower(x1) - slope * x1
            push!(pwl,LinearPiece(x1,x2,slope,b,ParaToFncLin(slope,b)))
            break
        end

        if f(x1) * f(x2) > 0 && abs(f(x1)) < EPS
            xp = perform_binary_sarch_lowf(f, x1, x2)
            # println("f(xp) = $(f(xp))")
            slope = du(xp)
            b = lower(x1) - slope * x1
            push!(pwl,LinearPiece(x1,xp,slope,b,ParaToFncLin(slope,b)))
            x1 = xp
            continue
        end

        a = find_zero(f, x1, x2)
        if isnan(a) && abs(f(x1)) < abs(f(x2)) a = x1
        elseif isnan(a) && abs(f(x2)) < abs(f(x1)) a = x2
        end
        
        @assert(!isnan(a), "a is still NAN!, f(x1) = $(f(x1)), f(x2) = $(f(x2))")

        slope = du(a)
        b = lower(x1) - slope * x1

        #find the second end of the segment
        if Δ(x2) <= 0.0
            # nextX1 = find_zero(Δ,(a,x2),Bisection(), atol = tol / 100)
            nextX1 = find_zero(Δ, a, x2)
        else
            nextX1 = x2
        end

        push!(pwl,LinearPiece(x1,nextX1,slope,b,ParaToFncLin(slope,b)))

        x1 = nextX1
    end
    # println("end of loop")
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