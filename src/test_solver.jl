using Plots
using CSV
using DataFrames
using JSON
include("ode_solver.jl")
include("evo_utils.jl")
include("network_cleanup.jl")
# include("settings.jl")

function gettopmodel(species_by_IDs)
    maxfitness = 0
    topnetwork = nothing
    topspecies = nothing
    for speciesID in keys(species_by_IDs)
        species = species_by_IDs[speciesID]
        if species.topfitness > maxfitness
            maxfitness = species.topfitness
            topnetwork = species.topnetwork
            topspecies = species
        end
    end
    #println("species: $(topspecies.ID)")
    #println("num offspring: $(topspecies.numoffspring)")
    return topnetwork, maxfitness
end


function main()
    pathtosettings = "/home/hellsbells/Desktop/networkEv/test_files/oscillator1.json"

    settings = read_usersettings(pathtosettings)

    objfunct = get_objectivefunction(settings)
    ng = get_networkgenerator(settings)

    # r = generate_random_reaction(ng)
    # print(r)
    # l = generate_reactionlist(ng)
    # println(l)
    # print(current_innovation_num_by_reaction)
    # population = evolve(settings, ng, objfunct)

    DELTA = .5
    NUM_GENERATION = 150

    population = generate_network_population(settings, ng)

    species_by_IDs = initialize_species_by_IDs(population)
    species_by_IDs = speciate(species_by_IDs, population, DELTA)

    species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs)
    species_by_IDs = calculate_num_offspring(species_by_IDs, total_fitness, settings)

    fitnesses = []
    # populationsizes = []
    for i in 1:NUM_GENERATION
        # fitness_by_species

        # global total_fitness
        # global numoffspring_by_species
        # global newpopulation
        if i%10 == 0
            println("starting generation $i")
        end
        println("starting generation $i")
        species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs)
        species_by_IDs = calculate_num_offspring(species_by_IDs, total_fitness, settings)

        bestnetwork, maxfitness = gettopmodel(species_by_IDs)

        push!(fitnesses, maxfitness)
        # println(maxfitness)
        
        

        population = reproduce_networks(species_by_IDs, settings, ng, objfunct, generation = i)
        
        # if i%10 == 0
        #     println("gen $i: num species $(length(keys(networks_by_species))) and pop size $(length(newpopulation))")
        #     bestnetwork, maxfitness = get_top_model(networks_by_species, objfunct)
        #     println("maxfitness: $maxfitness")
        #     astr =convert_to_antimony(bestnetwork)
        #     println(astr)

        # end
        # push!(populationsizes, length(newpopulation))


        species_by_IDs = speciate(species_by_IDs, population, DELTA)


        # println("This is generation $i and there are $(length(keys(networks_by_species)))")
    end
    # println(species_by_IDs)
    species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs)
    bestnetwork, maxfitness = gettopmodel(species_by_IDs)


    println("Best fitness: $maxfitness")
    astr = convert_to_antimony(bestnetwork)
    println(astr)
end


main()

# astr = "S0 -> S0 + S0; k1*S0
# S0 -> S2 + S2; k2*S0
# S0 + S2 -> S0; k3*S0*S2
# k1 = 146.49479609999202
# k2 = 37.73881143052887
# k3 = 132.0990695275749
# S0 = 1.0h
# S1 = 5.0
# S2 = 9.0"

# model = convert_from_antimony_string(astr)
# println(model.boundaryspecies)

# pathtosettings = "/home/hellsbells/Desktop/networkEv/test_files/oscillator1.json"

# settings = read_usersettings(pathtosettings)

# objfunct = get_objectivefunction(settings)

# sol = solve_ode(objfunct, model)

# print(sol.t)

# using Plots
# plot(sol)

# savefig("plot.png")

# println(populationsizes)

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