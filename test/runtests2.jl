using LinA
# f(x) = log(1 + (MathConstants.e - 1) * x)
# f(x) = 84.65625859999999/3*log(1 + (MathConstants.e^3 - 1)/392*x)
# pwl = LinA.Linearize(f, 0, 1, LinA.Relative(.0), LinA.ExactLin(); bounding = Under())
f(x) = 2.1525441506373877 * (sin(0.6666666666666666x) + 0.3333333333333333x) ^ 2
# pwl = LinA.Linearize(f, 0.0, 51.611965507001564, LinA.Relative(6.25), LinA.HeuristicLin(); bounding = Under())
pwl = Linearize(f, 0.45008006908541487, 15, Relative(50.0), ExactLin(); bounding = Under())
println("number of pieces = $(length(pwl))")
println.(pwl)