
"""
isContinuous(pwl, ε = 1e-5)

Determine whether a pwl function is continuous up to a numerical precision of ε.
# Arguments
- `plw` : pwl function
- `ε` : numerical precision

"""
function isContinuous(pwl, ε = 1e-5) 

    
    for i in 1:length(pwl)-1
        
        temp = pwl[i].xMax
        
        if pwl[i](temp) - pwl[i+1](temp) > ε
            return false
        end
        
    end
    
    true
    
end





"""
CplexBreakpoints(pwl, ε = 1e-5)

Outputs the breakpoints needed for CPLEX to natively model PWL functions.
# Arguments
- `plw` : pwl function
- `ε` : numerical precision

"""
function CplexBreakpoints(pwl, ε = 1e-5) 

    
    bpx = [pwl[1].xMin]

    
    for i in 1:length(pwl)-1
        
        temp = pwl[i].xMax
        push!(bpx,temp)
        
        if pwl[i](temp) - pwl[i+1](temp) > ε
            push!(bpx , temp)
        end
        
    end
    
    push!(bpx,pwl[end].xMax)
    
    return bpx
    
    
end