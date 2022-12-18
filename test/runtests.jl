using LinA
using Test

function testErrorAbsolute(f,pwl,delta)

    from = pwl[1].xMin
    to = pwl[end].xMax
    gap = (to-from)/1000

    diff(x) = abs(f(x) - pwl(x))

   all(diff.(from:gap:to) .<= delta + 1e-8)

end


@testset "Heuristic" begin

    pwl = Linearize(:(x^2),-10,10,Absolute(2))
    @test length(pwl) == 5
    @test testErrorAbsolute(x->x^2,pwl,2)

    @test length(Linearize(:(x^3),-10,10,Absolute(2),)) == 20


end

@testset "Exact" begin

    @test length(Linearize(:(x^2),-10,10,Absolute(2),ExactLin())) == 5
    @test length(Linearize(:(x^3),-10,10,Absolute(2),ExactLin())) == 19


end