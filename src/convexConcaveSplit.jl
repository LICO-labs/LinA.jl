function ConcavitySplit(x1::Real,x2::Real,expr_fnc::Ef)
    d2_f = Derive(Derive(expr_fnc))
    #special case for constant second derivative
    typeof(d2_f) <: Number && return Float64[]

    temp =  FctMaker(d2_f)
    f(x) = temp(x)

    try
        return find_zeros(f,x1,x2)
    catch y
        #if the second derivative is almost zero on an interval
        if isa(y, InexactError) || isa(y,ArgumentError)
            @warn "No concavity changes were found in the interval.\n It might be because of precision issues since the second derivative is almost constant"
            return Float64[]
        else
            rethrow()
        end
    end
end

