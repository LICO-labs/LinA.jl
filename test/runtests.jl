using LinA
using Test

function testErrorAbsolute(f,pwl,delta)

    from = pwl[1].xMin
    to = pwl[end].xMax
    gap = (to-from)/100

    diff(x) = abs(f(x) - pwl(x))

   all(diff.(from:gap:to) .<= delta + 1e-8)

end

function testErrorRelative(f,pwl,eps)

    from = pwl[1].xMin
    to = pwl[end].xMax
    gap = (to-from)/100

    diff(x) = abs(f(x) - pwl(x))

   all(diff.(from:gap:to) .<= eps * abs.(f.(from:gap:to)) .+ 1e-8)

end


@testset "Heuristic Absolute" begin

    pwl = Linearize(:(x^2),-10,10,Absolute(2))
    @test length(pwl) == 5
    @test testErrorAbsolute(x->x^2,pwl,2)

    pwl = Linearize(:(-x^2),-10,10,Absolute(2))
    @test length(pwl) == 5
    @test testErrorAbsolute(-x->x^2,pwl,2)

    pwl = Linearize(:(x^3),-10,10,Absolute(2))
    @test length(pwl) == 20
    @test  testErrorAbsolute(x->x^3,pwl,2)

    pwl = Linearize(x->x^3,-10,10,Absolute(2))
    @test length(pwl) == 20
    @test  testErrorAbsolute(x->x^3,pwl,2)

end

@testset "Heuristic Relative" begin

    pwl = Linearize(:(x^2+1),-10,10,Relative(1))
    @test length(pwl) == 22
    @test testErrorRelative(x->x^2+1,pwl,1)

    pwl = Linearize(:(-x^2-1),-10,10,Relative(1))
    @test length(pwl) == 22
    @test testErrorRelative(x->-x^2-1,pwl,1)

    pwl = Linearize(:(x^3+20),-2,2,Relative(1))
    @test length(pwl) == 6
    @test  testErrorRelative(x->x^3+20,pwl,1)

    pwl = Linearize(x->x^3+20,-2,2,Relative(1))
    @test length(pwl) == 6
    @test  testErrorRelative(x->x^3+20,pwl,1)

end




@testset "Exact Absolute" begin


    pwl = Linearize(:(x^2),-10,10,Absolute(2),ExactLin())
    @test length(pwl) == 5
    @test testErrorAbsolute(x->x^2,pwl,2)

    pwl = Linearize(:(-x^2),-10,10,Absolute(2),ExactLin())
    @test length(pwl) == 5
    @test testErrorAbsolute(x->-x^2,pwl,2)

    pwl = Linearize(:(x^3),-10,10,Absolute(2),ExactLin())
    @test length(pwl) == 19
    @test  testErrorAbsolute(x->x^3,pwl,2)

    pwl = Linearize(x->x^3,-10,10,Absolute(2),ExactLin())
    @test length(pwl) == 19
    @test  testErrorAbsolute(x->x^3,pwl,2)


end



@testset "Exact Relative" begin

    pwl = Linearize(:(x^2+1),-10,10,Relative(1),ExactLin())
    @test length(pwl) == 22
    @test testErrorRelative(x->x^2+1,pwl,1)

    pwl = Linearize(:(-x^2-1),-10,10,Relative(1),ExactLin())
    @test length(pwl) == 22
    @test testErrorRelative(x->-x^2-1,pwl,1)

    pwl = Linearize(:(x^3+20),-2,2,Relative(1),ExactLin())
    @test length(pwl) == 5
    @test  testErrorRelative(x->x^3+20,pwl,1)

    pwl = Linearize(x->x^3+20,-2,2,Relative(1),ExactLin())
    @test length(pwl) == 5
    @test  testErrorRelative(x->x^3+20,pwl,1)

end

@testset "Edge cases" begin

    #constant functions 

    pwl = Linearize(:(x-x),-10,10,Absolute(1))
    @test length(pwl) == 1

    
    pwl = Linearize(:(x-x),-10,10,Absolute(1),ExactLin())
    @test length(pwl) == 1

    #linear functions 

    pwl = Linearize(:(3x),-10,10,Absolute(1))
    @test length(pwl) == 1

    
    pwl = Linearize(:(3x),-10,10,Absolute(1),ExactLin())
    @test length(pwl) == 1


end

@testset "Exceptions" begin


    pwl = Linearize(:(x^2),-10,10,Absolute(1))
    @test_throws DomainError pwl(30)

    @test_throws ArgumentError Linearize(:(sin(x)^2),-Inf,2,Absolute(0.05))
end