

#module decoupage

#export decoupageConvavite

#using Roots
#using Calculus
#using GeneralizedGenerated


function seconde_derive(expr_fnc::Expr)
     d1_expr = differentiate(expr_fnc, :x)
     d1_expr = Calculus.simplify(d1_expr)
     d2_expr = differentiate(d1_expr, :x)
     d2_expr = Calculus.simplify(d2_expr)
    return d2_expr
end



 function decoupageConvavite(x1::Real,x2::Real,expr_fnc::Expr)

    d2_f = seconde_derive(expr_fnc)
    temp =  mk_function(:((x -> $d2_f)))
    f(x) = temp(x)
    #@eval f(x) = $d2_f

    try
        #return Base.invokelatest(find_zeros,f,x1,x2)
        return find_zeros(f,x1,x2)
    catch y
        #println(y)
        #si la dérivée est une constante find zeros lance une erreur
           if isa(y, InexactError)
            if differentiate(d2_f, :x) == 0
                return Float64[]
            end
            @warn "No concavity changes were found in the interval, might be because of precision issues"
           end
       end
end


#end
