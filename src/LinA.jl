 module LinA

using Roots
using Calculus
using GeneralizedGenerated



export Linearize , LinearBounding
export Absolute, Relative, Best, Under, Over, HeuristicLin, ExactLin

include("strucDef.jl")

include("Linearize Dispatch.jl")
#include("ORourke.jl")
#include("exact method.jl")

end # module
