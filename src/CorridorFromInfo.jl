




"""
    corridorFromInfo(x1::Real,x2::Real,expr_fct::Ef,e::ErrorType,bounding::BoundingType)

Create a corridor from an error type and bounding style. This corridor is a tuple (start,end,lowerBound, upperbound, derivative of the upper bound)

# Arguments
- `x1` : from
- `x2` : to
- `expr_fct` : function to linearize (either an expresion or a native julia function)
- `e` : error type. Both Absolute() and Relative() are implemented
- `bounding` : `Under()` for an underestimation, `Over()` for an overestimation, `Best()` for estimation that can go under or over the function.
"""
function corridorFromInfo(x1::Real,x2::Real,expr_fct::Ef,e::ErrorType,bounding::BoundingType)

    feval =  fctMaker(expr_fct)

    f(y::Real) = feval(y)

    temp = derive(expr_fct)
    tempEval = fctMaker(temp)
    

    df = y::Real -> tempEval(y)
    #note that upper return the function and it's derivative
    return x1,x2,lower(x1,x2,f,e,bounding),upper(x1,x2,f,df,e,bounding)...
end

#In the case of underestimation or overestimation one of the bound is the function itself
upper(x1::Real,x2::Real,f::Function,df::Function,e::ErrorType,::Under)  = x::Real -> f(x),df
lower(x1::Real,x2::Real,f::Function,e::ErrorType,::Over)  = x::Real -> f(x)


#in the case of best approximation, the bounds are defined the same as under for the lower function
#and upper for the top function
lower(x1::Real,x2::Real,f::Function,e::ErrorType,::Best)  = lower(x1,x2,f,e,Under())
upper(x1::Real,x2::Real,f::Function,df::Function,e::ErrorType,::Best) = upper(x1,x2,f,df,e,Over())


#absolute error case
lower(x1::Real,x2::Real,f::Function,e::Absolute,::Under)  = x::Real -> f(x) - e.delta
upper(x1::Real,x2::Real,f::Function,df::Function,e::Absolute,::Over) = x::Real -> f(x) + e.delta,df


#Relative error case
function lower(x1::Real,x2::Real,f::Function,e::Relative,::Under)
   # TODO ADD TOLERANCE as a parameter
   # This 10.0^(-5) seems arbitrairy but commes from the litterature

    if f((x1+x2)/2) >= 0
       return  x::Real -> f(x)*(1 - e.percent/100) - 10.0^(-5)
    end
    return x::Real -> f(x)*(1 + e.percent/100) - 10.0^(-5)
end

function upper(x1::Real,x2::Real,f::Function,df::Function,e::Relative,::Over)
   # TODO ADD TOLERANCE as a parameter
   # This 10.0^(-5) seems arbitrairy but commes from the litterature

    if f((x1+x2)/2) >= 0
       return x::Real -> f(x)*(1 + e.percent/100) + 10.0^(-5),x->(1 + e.percent/100)*df(x)
    end
    return x::Real -> f(x)*(1 - e.percent/100) + 10.0^(-5),x->(1 - e.percent/100)*df(x)
end
