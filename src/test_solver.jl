using Plots
using CSV
using DataFrames
using JSON
include("ode_solver.jl")
include("evo_utils.jl")
include("network_cleanup.jl")
# include("settings.jl")

using Dates
println("Starting $(now())")

using TimerOutputs

to = TimerOutput()

function gettopmodel(species_by_IDs::Dict{String, Species})
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
    pathtosettings = "/home/hellsbells/Desktop/networkEv/test_files/updownObjFunc.json"

    settings = read_usersettings(pathtosettings)

    objfunct = get_objectivefunction(settings)
    ng = get_networkgenerator(settings)

    DELTA = .65
    TARGET_NUM_SPECIES = 10
    SPECIES_MOD_STEP = 0.1
    NUM_GENERATION = 10#400
    break_gen = 100

    
    
    population = generate_network_population(settings, ng)

    
    species_by_IDs = initialize_species_by_IDs(population)
    species_by_IDs, DELTA = speciate(species_by_IDs, population, DELTA, TARGET_NUM_SPECIES, SPECIES_MOD_STEP)

    species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs)
    species_by_IDs = calculate_num_offspring(species_by_IDs, total_fitness, settings)

    fitnesses = []
    println("starting")

    for i in 1:NUM_GENERATION

        @timeit to "evaluate pop fitness" begin
            species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs)
        end
        @timeit to "calculate offspring" begin
            species_by_IDs = calculate_num_offspring(species_by_IDs, total_fitness, settings)
        end
        # L = length(keys(species_by_IDs))
        # println("$L species")
        # if L <=1
        #     println("generation $i: $L")
        # end
        

        bestnetwork, maxfitness = gettopmodel(species_by_IDs)
        # if maxfitness > 0.009
        #     println("generation $i and model $(bestnetwork.ID), $(maxfitness)")
        #     writeoutnetwork(bestnetwork, "model_$(bestnetwork.ID)_writeout", directory="final_models")

        # end

        push!(fitnesses, maxfitness)
    
        
        
        @timeit to "reproduce networks" begin
            population = reproduce_networks(species_by_IDs, settings, ng, objfunct, generation = i)
        end

        @timeit to "speciate" begin
            species_by_IDs, DELTA = speciate(species_by_IDs, population, DELTA, TARGET_NUM_SPECIES, SPECIES_MOD_STEP)
        end

        
    end
    
    species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs)
    bestnetwork, maxfitness = gettopmodel(species_by_IDs)


    println("Best fitness: $maxfitness")

    # for ID in keys(species_by_IDs)
    #     species = species_by_IDs[ID]
    #     println(species.topfitness)
    #     if species.topfitness > settings.writeout_threshold
    #         println("WRiting out runner up from species $(species.ID)")
    #         writeoutnetwork(species.topnetwork, "runnerup_$(species.ID)", directory="final_models")
    #     end
    # end

    
    astr = convert_to_antimony(bestnetwork)
    writeoutnetwork(bestnetwork, "model_$(bestnetwork.ID)", directory="final_models")
    # return fitnesses
end

main()
show(to)


# using Profile
# using ProfileView

# @profview main()

# @profview main()
# readlines()

# for i in 1:100
#     main()
    
# end

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