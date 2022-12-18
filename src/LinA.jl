 module LinA

using Roots
using Calculus
using GeneralizedGenerated



export Linearize , LinearBounding
export Absolute, Relative, Best, Under, Over, HeuristicLin, ExactLin
export isContinuous, CplexBreakpoints

include("strucDef.jl")
include("Utilities.jl")
include("convexConcaveSplit.jl")
include("CorridorFromInfo.jl")


include("ConvexCorridor.jl")
include("Heuristic.jl")

include("ORourke.jl")
include("exactPiece.jl")

include("exact method.jl")

include("Linearize Dispatch.jl")
include("Bounding.jl")
#include("ORourke.jl")
#include("exact method.jl")

end # module
