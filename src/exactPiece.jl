
include("ORourke.jl")

function exactPiece(start::Real,maximum::Real,lower,upper)
    
    #newMax = maximum
    line = LinearPiece(0,0,0,0,x->0)
    pts = collect(range(start,maximum,length=40))
    data = fctSample.(pts, lower,upper)
    
    succes=false;
    topIntersec = []
    lowIntersec = []
    
    while !succes
        
        crossing = true
        
        while crossing
            
            sort!(pts)
            data = fctSample.(pts, lower,upper)
            line = ORourke(data)
            
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
                if topDistance((topIntersec[i]+topIntersec[i+1])/2) < 0
                    push!(pts,(topIntersec[i]+topIntersec[i+1])/2)
                    crossing = true;
                end
                
            end
            
            for i in 1: length(lowIntersec) -1
                
                #other criteria if differentiable
                #if lowerDistance'(lowIntersec[i]) < 0
                
                if lowerDistance((lowIntersec[i] + lowIntersec[i+1])/2) < 0
                    push!(pts,(lowIntersec[i] + lowIntersec[i+1])/2)
                    push!(pts,lowIntersec[i])
                    crossing = true;
                end
                
            end
            
        end

        lastCovered = line.xMax
        
        if lastCovered == pts[end]
            return line
        end
        
        
        index = findfirst(isequal(lastCovered), pts)
        notCover = pts[index+1]
        
        
        if notCover - lastCovered < 0.0001 #|| notCover == newMax
            return line
        end
    
        
        pts = pts[1 : index + 1 ]
        
        
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
        #newMax = notCover
        
        
    end
    
 
end
