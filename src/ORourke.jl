
using Polyhedra, CDDLib

#2020/03/13
function pointPlane(p::dataError)
    # a and b are the variables to optimize
    # -x*a - b < - ymin
    # x*a + b < ymax
    return HalfSpace([-p.x, -1], -p.yMin) ∩ HalfSpace([p.x, 1], p.yMax)
end

function ORourke(pts)
    pts = copy(pts)
    
    poly = pointPlane(pts[1])
    temp = poly;
    coveredPts = [popfirst!(pts)];
    
   
    while !isempty(pts) && !isempty(polyhedron(poly, CDDLib.Library(:exact))) 
        
        temp = copy(poly);
        poly = poly ∩ pointPlane(pts[1])
        push!(coveredPts,popfirst!(pts)) 
    end
   
    
    
    if isempty(polyhedron(poly, CDDLib.Library(:exact)))
        pop!(coveredPts)
        poly = temp
    end

    
    #on prend un point quelconque dans le polygone de contraintes des coef
    verticesPoly = collect(points(doubledescription(poly)))

    lineCoef = sum(verticesPoly)/length(verticesPoly)
    
    #on le transforme en segement de droite
    line =  LinearPiece(coveredPts[1].x,coveredPts[end].x,lineCoef[1],lineCoef[2],x-> lineCoef[1]*x + lineCoef[2])
    
    return line
    
end