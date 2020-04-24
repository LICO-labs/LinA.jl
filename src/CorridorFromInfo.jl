
#include("strucDef.jl")
#using Calculus
#using GeneralizedGenerated

function corridorFromInfo(x1::Real,x2::Real,expr_fct::Expr,e::ErrorType,bounding::BoundingType)
   
    feval =  mk_function(:((x -> $expr_fct)))
   
    #@eval x::Real -> $expr_fct
    
    f(y::Real) = feval(y)
   
    temp = differentiate(expr_fct, :x)
    tempEval = mk_function(:((x -> $temp)))
    #tempEval = (@eval x::Real -> $temp)
    
    df = y::Real -> tempEval(y)
    #note that upper return the function and it's derivative
    return x1,x2,lower(x1,x2,f,e,bounding),upper(x1,x2,f,df,e,bounding)...
end

#In the case of underestimation or overestimation one of the bound is the function itself
upper(x1::Real,x2::Real,f::Function,df::Function,e::ErrorType,::Under)  = x::Real -> f(x),df
lower(x1::Real,x2::Real,f::Function,e::ErrorType,::Over)  = x::Real -> f(x)

#in the case of best approximation, the bound are defined the same as under for the lower function
#and upper for the top function
lower(x1::Real,x2::Real,f::Function,e::ErrorType,::Best)  = lower(x1,x2,f,e,Under()) 
upper(x1::Real,x2::Real,f::Function,df::Function,e::ErrorType,::Best) = upper(x1,x2,f,df,e,Over())

#absolute error case

lower(x1::Real,x2::Real,f::Function,e::Absolute,::Under)  = x::Real -> f(x) - e.delta
upper(x1::Real,x2::Real,f::Function,df::Function,e::Absolute,::Over) = x::Real -> f(x) + e.delta,df

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

