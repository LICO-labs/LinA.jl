

#module decoupage

#export decoupageConvavite

#using Roots
#using Calculus
#using GeneralizedGenerated


seconde_derive(fct::Ef) = derive(derive(fct))



 function decoupageConvavite(x1::Real,x2::Real,expr_fnc::Ef)

    d2_f = seconde_derive(expr_fnc)
    temp =  fctMaker(d2_f)
    f(x) = temp(x)
    #@eval f(x) = $d2_f
    #println(f(3),typeof(f))
    #return temp, f
    try
        #return Base.invokelatest(find_zeros,f,x1,x2)
        return find_zeros(f,x1,x2)
    catch y
        #println(y)
        #si la dérivée est une constante find zeros lance une erreur
           if isa(y, InexactError)
            #if differentiate(d2_f, :x) == 0
            #    return Float64[]
            #end
            @warn "No concavity changes were found in the interval.\n It might be because of precision issues since the second derivative is almost constant"
            return Float64[]
           end
       end
end


#end
