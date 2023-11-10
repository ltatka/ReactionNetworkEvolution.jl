using Plots
using CSV
using DataFrames
using JSON
include("ode_solver.jl")
include("evo_utils.jl")
include("network_cleanup.jl")
# include("settings.jl")

pathtosettings = "/home/hellsbells/Desktop/networkEv/src/test.json"

settings = read_usersettings(pathtosettings)

objfunct = get_objectivefunction(settings)
ng = get_networkgenerator(settings)

r = generate_random_reaction(ng)
print(r)
l = generate_reactionlist(ng)
println(l)
# print(current_innovation_num_by_reaction)
population = evolve(settings, ng, objfunct)
print("done")


# # ng2 = NetworkGenerator(settings.specieslist, settings.initialconditions, settings.nreactions,
# # settings.reactionprobabilities, settings.rateconstantrange, seed=population[1])
# # population = evolve(settings, ng, objfunct)

# bestnetwork = population[1]
# # println(convert_to_antimony(bestnetwork))
# bestnetwork = cleanupreactions(bestnetwork)
# println(convert_to_antimony(bestnetwork))

# solution = solve_ode(objfunct, bestnetwork)

# # # display(solution.u)

# using Plots
# plt = plot(solution)
# savefig(plt, "/home/hellsbells/Desktop/attemptedoscillator2.png")
# println("SUCCESS")
# println(NUM_VARIANTS)
# plt2 = plot(NUM_VARIANTS)
# savefig(plt2, "/home/hellsbells/Desktop/variants.png")