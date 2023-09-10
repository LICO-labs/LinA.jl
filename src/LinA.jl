 module LinA

# using Roots
using Calculus
using GeneralizedGenerated
using IntervalArithmetic, IntervalRootFinding


export Linearize , LinearBounding
export Absolute, Relative, Best, Under, Over, HeuristicLin, ExactLin
export isContinuous, breakpoints

include("strucDef.jl")
include("utilities.jl")
include("convexConcaveSplit.jl")
include("corridorFromInfo.jl")


include("convexCorridor.jl")
include("heuristic.jl")

include("oRourke.jl")
include("exactPiece.jl")

include("exactMethod.jl")

include("linearizeDispatch.jl")
include("bounding.jl")


end # module
