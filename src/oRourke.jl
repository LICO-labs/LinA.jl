


function PointPlane(p::dataError)
    # a and b are the variables to optimize
    # -x*a - b < - ymin
    # x*a + b < ymax
    return HalfSpace([-p.x, -1], -p.yMin) ∩ HalfSpace([p.x, 1], p.yMax)
end

function ORourke(pts)
    pts = copy(pts)

    poly = PointPlane(pts[1])
    temp = poly
    coveredPts = [popfirst!(pts)]

    #add points until no line is feasible or all pts are covered
    while !isempty(pts) && !isempty(polyhedron(poly, CDDLib.Library(:exact)))
        temp = copy(poly)
        poly = poly ∩ PointPlane(pts[1])
        push!(coveredPts, popfirst!(pts))

    end



    if isempty(polyhedron(poly, CDDLib.Library(:exact)))
        pop!(coveredPts)
        poly = temp
    end


    #We take any arbitrary point in the polyhedra (here the "center")
    verticesPoly = collect(points(doubledescription(poly)))
    lineCoef = sum(verticesPoly) / length(verticesPoly)

    line = LinearPiece(
        coveredPts[1].x,
        coveredPts[end].x,
        lineCoef[1],
        lineCoef[2],
        x -> lineCoef[1] * x + lineCoef[2],
    )

    return line

end
