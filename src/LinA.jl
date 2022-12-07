 module LinA

using Roots
using Calculus
using GeneralizedGenerated



export Linearize , LinearBounding
export Absolute, Relative, Best, Under, Over, HeuristicLin, ExactLin
export isContinuous, CplexBreakpoints

include("strucDef.jl")
include("Bounding.jl")
include("Linearize Dispatch.jl")
include("Utilities.jl")
#include("ORourke.jl")
#include("exact method.jl")

end # module
