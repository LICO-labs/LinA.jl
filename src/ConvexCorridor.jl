#using Roots
#include("strucDef.jl")

"""
    LinearizeConvex(x1,x2,lower::Function,upper::Function,du::Function)

Makes an optimal piecewise Linear approximation from x1 to x2 of a convex corridor
# Arguments
- lower : function at the bottom of the corridor
- upper : function on top of the corridor
- du : derivative of the upper function

"""
function LinearizeConvex(x1,x2,lower::Function,upper::Function,du::Function)

    #TODO mettre la tolérence en paramètre
    tol = 0.00001 #10.0^(-5)
    slope = 0.0;
    b = 0.0;
    pwl = Array{LinearPiece}(undef, 0)
    #x1 = round(x1_initial, digits = digits)

    lin(x) = slope * x + b
    #Δ represent the distance between the linear function and the bottom of the corridor
    Δ(x) = lin(x) - lower(x)
    #distance between the linear function and the bottom of the corridor at the start of the coridor if tangeant in "a"
    f(a)= upper(a) + du(a)*(x1-a) - lower(x1)

    while x2 - x1 > tol

        #find the tangence point if possible
        if f(x2) >= 0
            slope = (upper(x2) - lower(x1)) / (x2 - x1)
            b = lower(x1) - slope * x1
            push!(pwl,LinearPiece(x1,x2,slope,b,paraToFncLin(slope,b)))
            break
        end

        a = find_zero(f,(x1,x2), Bisection())
        slope = du(a)
        b = lower(x1) - slope * x1

        #find the second end of the segment
        if Δ(x2) <= 0
           nextX1 = find_zero(Δ,(a,x2),Bisection())
        else
            nextX1 = x2
        end

        push!(pwl,LinearPiece(x1,nextX1,slope,b,paraToFncLin(slope,b)))

        x1 = nextX1
    end

    return pwl

end
