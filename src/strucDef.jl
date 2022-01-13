

import LinearAlgebra.adjoint , Base.-, Base.+, Base.show

using Calculus, ForwardDiff

const Ef=Union{Expr, Function}

#to use the ' notattion for derivative (which was standard notation in Calculus.jl)
adjoint(f::Function) = derivative(f)

-(expr::Expr) = Expr(:call, :-, expr)
-(f::Function) = x->-f(x)
#-(expr::Expr) = :(($expr) * -1)
paraToFncLin(a,b) = x -> a*x + b


struct LinearPiece <: Function
    xMin::Real
    xMax::Real
    a::Real
    b::Real
    fct::Function
end


function Base.show(io::IO,::MIME"text/plain", m::LinearPiece)
   print(io,m.a," x ",m.b >= 0 ? "+ " : "",m.b," from ",m.xMin," to ",m.xMax)
end


(l::LinearPiece)(x::Real) = l.fct(x)
-(p::LinearPiece) = LinearPiece(p.xMin,p.xMax,-p.a,-p.b,paraToFncLin(-p.a,-p.b))
+(p::LinearPiece,x::Real) = LinearPiece(p.xMin,p.xMax,p.a,p.b + x,paraToFncLin(p.a,p.b + x))
-(p::LinearPiece,x::Real) = LinearPiece(p.xMin,p.xMax,p.a,p.b - x,paraToFncLin(p.a,p.b - x))



#evaluate a pieceWise linear function

function (pwl::Array{LinearPiece, 1})(x::Real)

    if x < pwl[1].xMin || x > pwl[end].xMax
        throw(DomainError(x, "argument must be in the domain of the function"))
    end

    starts = getfield.(pwl,:xMin)
    pieceIndex = searchsortedlast(starts,x)
    return pwl[pieceIndex](x)

end

+(pwl::Array{LinearPiece, 1}, x::Real) = pwl .+ x
-(pwl::Array{LinearPiece, 1}, x::Real) = pwl .- x


abstract type ErrorType end

struct Absolute <: ErrorType
    delta::Real
end

struct Relative <: ErrorType
    percent::Real
end

abstract type BoundingType end

struct Best <:BoundingType end
struct Under <:BoundingType end
struct Over <:BoundingType end
struct UnderOver <:BoundingType end

-(::Under) = Over()
-(::Over) = Under()
-(x::Best) = x


abstract type AbstractAlgorithm end
struct HeuristicLin <:AbstractAlgorithm end
struct ExactLin <:AbstractAlgorithm end

#For the exact method

struct dataError
    x::Real
    yMin::Real
    yMax::Real
end

function fctSample(x, lower, upper) 

    return dataError(x,lower(x),upper(x))
end

#for polymorphism between symbolic and not symbolic
fctMaker(e::Expr) = mk_function(:((x -> $e)))
fctMaker(e::Function) = e
fctMaker(x::Number) = y->x
derive(x::Number) = 0
derive(expr::Expr) = Calculus.simplify(differentiate(expr, :x))
derive(f::Function) = z->ForwardDiff.gradient(x->f(x[1]),[z])[1]

;
