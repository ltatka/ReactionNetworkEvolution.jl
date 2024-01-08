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

# r = generate_random_reaction(ng)
# print(r)
# l = generate_reactionlist(ng)
# println(l)
# print(current_innovation_num_by_reaction)
# population = evolve(settings, ng, objfunct)

DELTA = 0.2

population = generate_network_population(settings, ng)
networks_by_species = Dict(population[1].ID => population)
println(keys(networks_by_species))

networks_by_species, fitness_by_species, total_fitness = evaluate_population_fitness(objfunct, networks_by_species)
numoffspring_by_species = calculate_num_offspring(fitness_by_species, total_fitness, settings)


newpopulation = reproduce_networks(networks_by_species, numoffspring_by_species, settings, ng, objfunct)

new_networks_by_species = speciate(networks_by_species, newpopulation, DELTA)
println(keys(new_networks_by_species))

# newpopulation = []
# for network in population
#     network = mutatenetwork(settings, ng, network)
#     push!(newpopulation, network)
# end

# networks_by_species = speciate(networks_by_species, newpopulation, 0.2)
# println(keys(networks_by_species))

# for species in keys(networks_by_species)
#     a = sort(networks_by_species, by=fitness)

# population2 = generate_network_population(settings, ng)
# print(population2[1])

# network1 = population[1]
# network2 = population2[1]

# for r in network1.reactionlist
#     println(r)
# end
# println("*********************")

# for r in network2.reactionlist
#     println(r)
# end

# newreactions = crossover(network1, network2)
# println("******")
# for r in newreactions
#     println(r)
# end


# deltas = [0.1, .15: 1/2 ]

# for delta in deltas
#     new_networks_by_species = speciate(networks_by_species, newpopulation, delta)
#     println("$delta: $(length(new_networks_by_species))")

# end


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


print("done")