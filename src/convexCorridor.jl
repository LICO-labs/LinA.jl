

"""
    LinearizeConvex(x1,x2,lower::Function,upper::Function,du::Function)

Makes an optimal piecewise Linear approximation from x1 to x2 of a convex corridor
# Arguments
- lower : function at the bottom of the corridor
- upper : function on top of the corridor
- du : derivative of the upper function

"""
function LinearizeConvex(x1,x2,lower::Function,upper::Function,du::Function)

    
    tol = 1e-7
    slope = 0.0;
    b = 0.0;
    pwl = Array{LinearPiece}(undef, 0)
    

    lin(x) = slope * x + b
    #Δ represent the distance between the linear function and the bottom of the corridor
    Δ(x) = lin(x) - lower(x)
    relΔ(x) = Δ(x) / max(1e-2, x, -x)
    #distance between the linear function and the bottom of the corridor at the start of the coridor if tangeant in "a"
    f(a) = upper(a) + du(a)*(x1-a) - lower(x1)
    relf(a) = f(a) / max(1e-2, a, -a)

    # println(stderr, "loop")
    while x2 - x1 > tol
        # println(stderr, "asdasd")
        # println("entering loop with x1 = $x1, x2 = $x2, f(x1) = $(f(x1)), f(x2) = $(f(x2))")
        # if f(x2) < 0.0 && f(x1) * f(x2) > 0.0 break end
        #find the tangence point if possible
        if relf(x2) >= 0.0
            slope = (upper(x2) - lower(x1)) / (x2 - x1)
            b = lower(x1) - slope * x1
            push!(pwl,LinearPiece(x1,x2,slope,b,ParaToFncLin(slope,b)))
            break
        end

        # if f(x2) < 0.0 && f(x1) * f(x2) > 0.0 break end
        # println("finding zero with f(x1) = $(f(x1)), f(x2) = $(f(x2))")
        # if abs(f(x1)) < tol a = x1
        # elseif abs(f(x2)) < tol a = x2
        # else a = find_zero(f,(x1,x2), Bisection())
        # end
        # a = find_zero(f,(x1,x2), Bisection(), atol = tol/ 100)
        a = find_zero(relf, x1, x2)
        if isnan(a) && abs(f(x1)) < 1e-5 a = x1
        elseif isnan(a) && abs(f(x2)) < 1e-5 a = x2
        end
        
        @assert(!isnan(a), "a is still NAN!, relf(x1) = $(relf(x1)), relf(x2) = $(relf(x2)), f(x1) = $(f(x1)), f(x2) = $(f(x2))")

        slope = du(a)
        b = lower(x1) - slope * x1

        #find the second end of the segment
        if relΔ(x2) <= 0.0
            # nextX1 = find_zero(Δ,(a,x2),Bisection(), atol = tol / 100)
            nextX1 = find_zero(relΔ, a, x2)
        else
            nextX1 = x2
        end

        push!(pwl,LinearPiece(x1,nextX1,slope,b,ParaToFncLin(slope,b)))

        x1 = nextX1
    end

    return pwl

end
