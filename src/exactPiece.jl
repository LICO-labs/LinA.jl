


"""
    ExactPiece(start::Real,maximum::Real,lower,upper)

Computes the maximal linear piece starting at `start` which lies in between `lower` and `upper`. Works for any continuous `lower` and `upper`.
# Arguments
- `start` : from
- `maximum` : maximal end point of the linear segment
- `lower` : lower bound of the corridor
- `upper` : upper bound of the corridor
"""
function ExactPiece(start::Real,maximum::Real,lower,upper)
    #TODO: add epsilon as an argument for the user
    #TODO: If intersections are epsilon close skip intersections
    
    #numerical precision 
    epsilon = 1e-5 
    line = LinearPiece(0,0,0,0,x->0)
    pts = collect(range(start,maximum,length=50))
    data = FctSample.(pts, lower,upper)
    
    succes=false;
    topIntersec = []
    lowIntersec = []
    
    while !succes
        
        crossing = true
        
        while crossing
            
            sort!(pts)
            data = FctSample.(pts, lower,upper)
            line = ORourke(data)

            #find if the solution on the discretized problem works on the original problem 
            topDistance = x-> upper(x) - line.fct(x)
            lowerDistance = x-> line.fct(x) - lower(x)
            
            #try catch to handle rare cases with function asymptotic to zero
            try
                topIntersec = find_zeros(topDistance,line.xMin,line.xMax)
            catch
                topIntersec = []
            end
            
            try
                lowIntersec = find_zeros(lowerDistance,line.xMin,line.xMax)
            catch
                lowIntersec = []
            end

            
            crossing = false
            
            for i in 1: length(topIntersec) -1
                
                #other criteria if differentiable
                #if topDistance'(topIntersec[i]) < 0
                if topDistance((topIntersec[i]+topIntersec[i+1])/2) <- epsilon
                    
                    push!(pts,topIntersec[i])
                    push!(pts,(topIntersec[i]+topIntersec[i+1])/2)
                    push!(pts,topIntersec[i+1])
                    
                    #previously any precision of 1e-5 or below very rarely caused an infinite loop here because of the conversion 
                    #Rational{BigInt} <->  float64 used in ORourke ( method CDDLib.Library(:exact)) which is why randomization was used
                    push!(pts,RandomMidPoint(topIntersec[i], topIntersec[i+1]))
                    
                    crossing = true;
                end
                
            end
            
            for i in 1: length(lowIntersec) -1
                #other criteria if differentiable
                #if lowerDistance'(lowIntersec[i]) < 0
                if lowerDistance((lowIntersec[i] + lowIntersec[i+1])/2) < - epsilon
                    
                    
                    push!(pts,lowIntersec[i])
                    push!(pts,(lowIntersec[i] + lowIntersec[i+1])/2)
                    push!(pts,lowIntersec[i+1])
                
                end
                
            end
            
            
        end


        lastCovered = line.xMax
        
        if lastCovered == pts[end]
            return line
        end
        
        
        index = findfirst(isequal(lastCovered), pts)
        notCover = pts[index+1]

        #verify if 
        if notCover - lastCovered < epsilon #|| notCover == newMax
            return line
        end
    
        
        pts = pts[1 : index + 1 ]
        
        # Heuristic to achieve faster convergence (try exenting the segment until goes out of corridor)
        lExtend = maximum
        uExtend = maximum
        
        try 
            lExtend = find_zeros(x-> line.fct(x) - lower(x), line.xMax,maximum)[1]
            catch y
        end
        try 
            uExtend = find_zeros(x-> line.fct(x) - upper(x), line.xMax,maximum)[1]
            catch y
        end
        furthest = min(uExtend,lExtend)
        push!(pts, furthest)
        push!(pts, (notCover + furthest)/2 )
        push!(pts, (notCover + lastCovered)/2 )

        
        
    end
    
 
end
