

import LinearAlgebra.adjoint , Base.-, Base.+, Base.show

using Calculus, ForwardDiff

const Ef=Union{Expr, Function}


Minus(expr::Expr) = Expr(:call, :-, expr)
Minus(f::Function) =  x->-f(x)

ParaToFncLin(a,b) = x -> a*x + b


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
-(p::LinearPiece) = LinearPiece(p.xMin,p.xMax,-p.a,-p.b,ParaToFncLin(-p.a,-p.b))
+(p::LinearPiece,x::Real) = LinearPiece(p.xMin,p.xMax,p.a,p.b + x,ParaToFncLin(p.a,p.b + x))
-(p::LinearPiece,x::Real) = LinearPiece(p.xMin,p.xMax,p.a,p.b - x,ParaToFncLin(p.a,p.b - x))



#evaluate a pieceWise linear function

function (pwl::Array{LinearPiece, 1})(x::Real, tiebreak = Under())

    eps = EPS
    if x < pwl[1].xMin - eps || x > pwl[end].xMax + eps
        throw(DomainError(x, "argument must be in the domain of the function"))
    end

    starts = getfield.(pwl, :xMin)
    ends = getfield.(pwl, :xMax)
    piece_start = max(1, searchsortedlast(starts, x - eps))
    piece_end = min(length(pwl), searchsortedfirst(ends, x + eps))
    # println("piece_start = $piece_start, piece_end = $piece_end")
    if tiebreak isa Under
        return minimum([pwl[i](x) for i in piece_start:piece_end if x >= pwl[i].xMin - eps && x <= pwl[i].xMax + eps])
    elseif tiebreak isa Over
        return maximum([pwl[i](x) for i in piece_start:piece_end if x >= pwl[i].xMin - eps && x <= pwl[i].xMax + eps])
    else 
        midpiece = round(Int64, (piece_start + piece_end) / 2) 
        return pwl[midpiece](x)
    end

end

+(pwl::Array{LinearPiece, 1}, x::Real) = pwl .+ x
-(pwl::Array{LinearPiece, 1}, x::Real) = pwl .- x


abstract type ErrorType end

struct Absolute <: ErrorType
    delta::Real
end

"""
Relative error from the function (in percentage)

    !!! warning
    For a relative error to be well defined, the function needs to have no zeros!

"""
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

function FctSample(x, lower, upper) 

    return dataError(x,lower(x),upper(x))
end

function RandomMidPoint(x::Number,y::Number)
    r = rand()
    r * x + (1-r) * y
end

#for polymorphism between symbolic and not symbolic
FctMaker(e::Expr) = mk_function(:((x -> $e)))
FctMaker(e::Function) = e
FctMaker(x::Number) = y->x
Derive(x::Number) = 0
Derive(expr::Expr) = Calculus.simplify(differentiate(expr, :x))
Derive(f::Function) = z->ForwardDiff.gradient(x->f(x[1]),[z])[1]
#ForwardDiff.derivative(f, x::Real)

const EPS = 1e-9
