




"""
    CorridorFromInfo(x1::Real,x2::Real,expr_fct::Ef,e::ErrorType,bounding::BoundingType)

Create a corridor from an error type and bounding style. This corridor is a tuple (start,end,lowerBound, upperbound, derivative of the Upper bound)

# Arguments
- `x1` : from
- `x2` : to
- `expr_fct` : function to linearize (either an expresion or a native julia function)
- `e` : error type. Both Absolute() and Relative() are implemented
- `bounding` : `Under()` for an underestimation, `Over()` for an overestimation, `Best()` for estimation that can go under or over the function.
"""
function CorridorFromInfo(
    x1::Real,
    x2::Real,
    expr_fct::Ef,
    e::ErrorType,
    bounding::BoundingType,
)

    feval = FctMaker(expr_fct)

    f(y::Real) = feval(y)

    temp = Derive(expr_fct)
    TempEval = FctMaker(temp)


    df(y::Real) = TempEval(y)
    #note that Upper return the function and it's derivative
    return x1, x2, Lower(x1, x2, f, e, bounding), Upper(x1, x2, f, df, e, bounding)...
end

#In the case of underestimation or overestimation one of the bound is the function itself
Upper(x1::Real, x2::Real, f::Function, df::Function, e::ErrorType, ::Under) = f, df
Lower(x1::Real, x2::Real, f::Function, e::ErrorType, ::Over) = f

struct Shift{F,T}
    f::F
    Δ::T
end

(s::Shift)(x) = s.f(x) + s.Δ

#in the case of best approximation, the bounds are defined the same as under for the Lower function
#and Upper for the top function
Lower(x1::Real, x2::Real, f::Function, e::ErrorType, ::Best) = Lower(x1, x2, f, e, Under())
Upper(x1::Real, x2::Real, f::Function, df::Function, e::ErrorType, ::Best) =
    Upper(x1, x2, f, df, e, Over())

#absolute error case
Lower(x1::Real, x2::Real, f::Function, e::Absolute, ::Under) = Shift(f, -e.delta)
Upper(x1::Real, x2::Real, f::Function, df::Function, e::Absolute, ::Over) =
    Shift(f, +e.delta), df

struct Scale{F,T}
    f::F
    α::T
end

(s::Scale)(x) = s.f(x) * s.α

#Relative error case
function Lower(x1::Real, x2::Real, f::Function, e::Relative, ::Under; ε = 1e-5)
    # TODO ADD TOLERANCE as a parameter to the user
    # This 10.0^(-5) seems arbitrairy but commes from the litterature 
    # and it helps to gives a sensible definition for a relative corridor that starts at y=0 

    if f((x1 + x2) / 2) >= 0
        return Shift(Scale(f, 1 - e.percent / 100), -ε)
    end
    return Shift(Scale(f, 1 + e.percent / 100), +ε)
end

function Upper(x1::Real, x2::Real, f::Function, df::Function, e::Relative, ::Over; ε = 1e-5)
    # TODO ADD TOLERANCE as a parameter to the user
    # This 10.0^(-5) seems arbitrairy but commes from the litterature 
    # and it helps to gives a sensible definition for a relative corridor that starts at y=0 

    if f((x1 + x2) / 2) >= 0
        return Shift(Scale(f, 1 + e.percent / 100), +ε), Scale(df, 1 + e.percent / 100)#f(x)*(1 + e.percent/100) + 10.0^(-5),x->(1 + e.percent/100)*df(x)
    end
    return Shift(Scale(f, 1 - e.percent / 100), -ε), Scale(df, 1 - e.percent / 100)#f(x)*(1 - e.percent/100) + 10.0^(-5),x->(1 - e.percent/100)*df(x)
end
