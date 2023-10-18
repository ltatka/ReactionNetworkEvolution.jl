using Plots
using CSV
using DataFrames
using JSON
include("ode_solver.jl")
include("evo_utils.jl")
# include("settings.jl")

pathtosettings = "/home/hellsbells/Desktop/networkEv/src/test.json"

settings = read_usersettings(pathtosettings)

objfunct = get_objectivefunction(settings)
ng = get_networkgenerator(settings)

population = evolve(settings, ng, objfunct)

bestnetwork = population[1]
println(convert_to_antimony(bestnetwork))
bestnetwork = cleanupreactions(bestnetwork)
println(convert_to_antimony(bestnetwork))

solution = solve_ode(objfunct, bestnetwork)

display(solution.u)

using Plots
plt = plot(solution)
savefig(plt, "/home/hellsbells/Desktop/attemptedoscillator.png")

println("Success")
