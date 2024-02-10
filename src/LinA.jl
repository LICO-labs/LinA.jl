 module LinA

# using Roots
using Calculus
using GeneralizedGenerated
using IntervalArithmetic, IntervalRootFinding
using PrecompileTools
using Optim


export Linearize , LinearBounding, SimultaneousLin
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
include("simultaneousLin.jl")

@setup_workload begin
    f = x -> x * (x + 1) * log(x + 2) + 1
    @compile_workload begin
        for e in [Relative(1e-1), Absolute(1e-1)], alg in [HeuristicLin, ExactLin], bounding in [Under, Over, Best]
            pwl = Linearize(f, 0, 1, e, alg(); bounding = bounding())
            pwl(0.5, bounding)
        end
    end
end
end # module
