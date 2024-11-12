


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
    line = LinearPiece(0,0,0,0,x->0)
    pts = collect(range(start,maximum,length=50))
    data = FctSample.(pts, lower,upper)
    
    success = false
    topIntersec = []
    lowIntersec = []
    
    while !success
        
        crossing = true
        
        while crossing
            
            sort!(pts)
            data = FctSample.(pts, lower,upper)
            line = ORourke(data)

            #find if the solution on the discretized problem works on the original problem 
            topDistance = x-> upper(x) - line.fct(x)
            lowerDistance = x-> line.fct(x) - lower(x)

            topDistanceRel = x -> (topDistance(x) / max(EPS, abs(upper(x))))
            lowerDistanceRel = x -> (lowerDistance(x) / max(EPS, abs(lower(x))))
    
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
                midpoints = collect(topIntersec[i] : dp / 20 : topIntersec[i + 1])
                topdpoints = [(topDistanceRel(p), p) for p in midpoints]
                sort!(topdpoints)
                for (topd, p) in topdpoints[1 : 1]
                    if topd < - EPS
                        push!(pts,p)
                        crossing = true;
                    end
                end
            end
            
            for i in 1: length(lowIntersec) -1
                dp = lowIntersec[i + 1] - lowIntersec[i]
                midpoints = collect(lowIntersec[i] : dp / 20 : lowIntersec[i + 1])
                lowdpoints = [(lowerDistanceRel(p), p) for p in midpoints]
                sort!(lowdpoints)
                for (lowd, p) in lowdpoints[1 : 1]
                    if lowd < - EPS
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
        # index = findfirst(t -> abs(t - lastCovered) < 1e-7, pts)
        notCover = pts[index+1]
        @assert(!isnan(notCover), "notCover is NaN, points are $pts and index + 1 = $(index + 1)")

        #verify if 
        if notCover - lastCovered < EPS #|| notCover == newMax
            # println("breaking loop, notCover = $notCover, lastCovered = $lastCovered")
            return line
        else
            # println("continuing loop, notCover = $notCover, lastCovered = $lastCovered, $(notCover - lastCovered)")
        end
    
        
        pts = pts[1 : index + 1 ]
        
        # Heuristic to achieve faster convergence (try exenting the segment until goes out of corridor)
        lExtend = maximum
        uExtend = maximum
        
        try 
            lowerDistance = x -> line.fct(x) - lower(x)
            lowerDistanceRel = x -> lowerDistance(x) / max(EPS, abs(lower(x)))
            lExtend = find_zeros(lowerDistanceRel, line.xMax, maximum)[1]
            catch y
        end
        try 
            topDistance = x-> upper(x) - line.fct(x)
            topDistanceRel = x -> topDistance(x) / max(EPS, abs(upper(x)))
            uExtend = find_zeros(topDistanceRel, line.xMax, maximum)[1]
            catch y
        end
        @assert(!isnan(uExtend), "uExtend is nan")
        @assert(!isnan(lExtend), "lExtend is nan")
        @assert(!isnan(notCover), "notCover is nan")
        furthest = min(uExtend,lExtend)
        push!(pts, furthest)
        push!(pts, (notCover + furthest) / 2 )
        push!(pts, (notCover + lastCovered) / 2 )

        
        
    end
    
 
end
