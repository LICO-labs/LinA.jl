


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
    epsilon = 1e-7
    line = LinearPiece(0,0,0,0,x->0)
    pts = collect(range(start,maximum,length=50))
    data = FctSample.(pts, lower,upper)
    
    succes=false;
    topIntersec = []
    lowIntersec = []
    
    #find if the solution on the discretized problem works on the original problem 
    topDistance = x-> upper(x) - line.fct(x)
    lowerDistance = x-> line.fct(x) - lower(x)
    
    topDistanceRel = x-> topDistance(x) / max(1e-10, abs(upper(x)) + abs(line.fct(x)))
    lowerDistanceRel = x-> lowerDistance(x) / max(1e-10, abs(lower(x)) + abs(line.fct(x)))

    while !succes
        
        crossing = true
        
        while crossing
            
            sort!(pts)
            data = FctSample.(pts, lower,upper)
            line = ORourke(data)

            #try catch to handle rare cases with function asymptotic to zero
            try
                topIntersec = find_zeros(topDistanceRel,line.xMin,line.xMax)
            catch
                topIntersec = []
            end
            
            try
                lowIntersec = find_zeros(lowerDistanceRel,line.xMin,line.xMax)
            catch
                lowIntersec = []
            end

            
            crossing = false
            
            for i in 1: length(topIntersec) -1
                dp = topIntersec[i + 1] - topIntersec[i]
                midpoints = collect(topIntersec[i] : dp / 10 : topIntersec[i + 1])
                topdpoints = [(topDistanceRel(p), p) for p in midpoints]
                sort!(topdpoints)
                # n = ceil(Int64, length(topdpoints) / 5)
                # topdpoints = topdpoints[1 : n]
                for (topd, p) in topdpoints[1 : 1]
                    #other criteria if differentiable
                    #if topDistance'(topIntersec[i]) < 0
                    if topd < - epsilon
                        # println("adding point $p with violation of $topd)")
                        push!(pts,p)
                        crossing = true;
                    end
                end
            end
            
            for i in 1: length(lowIntersec) -1
                dp = lowIntersec[i + 1] - lowIntersec[i]
                midpoints = collect(lowIntersec[i] : dp / 10 : lowIntersec[i + 1])
                lowdpoints = [(lowerDistanceRel(p), p) for p in midpoints]
                sort!(lowdpoints)
                # n = ceil(Int64, length(lowdpoints) / 5)
                # lowdpoints = lowdpoints[1 : n]
                for (lowd, p) in lowdpoints[1 : 1]
                    #other criteria if differentiable
                    #if lowerDistance'(lowIntersec[i]) < 0
                    if lowd < - epsilon
                        push!(pts,p)
                        crossing = true
                    end
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
            lExtend = find_zero(lowerDistanceRel, line.xMax,maximum)
            catch y
        end
        try 
            uExtend = find_zero(topDistanceRel, line.xMax,maximum)
            catch y
        end
        furthest = min(uExtend,lExtend)
        push!(pts, furthest)
        push!(pts, (notCover + furthest)/2 )
        push!(pts, (notCover + lastCovered)/2 )

        
        
    end
    
 
end
