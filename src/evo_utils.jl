# include("reaction_network.jl")
include("settings.jl")
using Distributions
using CSV
using DataFrames






function addreaction(ng::NetworkGenerator, network::ReactionNetwork)
    reaction = generate_random_reaction(ng)
    push!(network.reactionlist, reaction)
    return network
end

function deletereaction(network::ReactionNetwork)
    idx = sample(1:length(network.reactionlist), 1)
    deleteat!(network.reactionlist, idx)
    return network
end

function adddeletereaction(ng:: NetworkGenerator, network::ReactionNetwork)
    p = rand()
    if p < 0.5
        network = addreaction(ng, network)
    else
        network = deletereaction(network)
    end
    return network
end


function mutaterateconstant(network::ReactionNetwork)
    idx = sample(1:length(network.reactionlist), 1)
    percentchange = rand(Uniform(.8, 1.2))
    network.reactionlist[idx][1].rateconstant *= percentchange
    return network
end

function mutatenetwork(settings::Settings, ng::NetworkGenerator, network::ReactionNetwork)
    p = rand()
    if p < settings.mutationprobabilities.adddeletereaction
        network = adddeletereaction(ng, network)
    else
        mutaterateconstant(network)
    end
    return network
end

function mutate_nonelite_population(settings::Settings, ng::NetworkGenerator, population)
    nelite = Int(floor(settings.portionelite*length(population)))
    for i in (nelite+1):length(population)
        population[i] = mutatenetwork(settings, ng, population[i])
    end
    return population
end

function generate_network_population(settings::Settings, ng::NetworkGenerator)
    population = []
    for i in 1:settings.usersettings.populationsize
        network = generate_random_network(ng)
        push!(population, network)
    end
    return population
end

function sortbyfitness(population)
    sort!(population, by=get_fitness)
    return population
end


function eliteselect(settings::Settings, population)
    nelite = Int(floor(settings.portionelite*length(population)))
    newpopulation = []
    for i in 1:nelite
        push!(newpopulation, deepcopy(population[i]))
    end
    return newpopulation
end


function tournamentselect(settings::Settings, population, newpopulation)
    #For now, this is only going to select networks and put them in the new population. Will mutate them later
    nelite = Int(floor(settings.portionelite*length(population)))
    for i in 1:(settings.usersettings.populationsize - nelite)
        # This will allow elites to be selected and they might dominate
        # What if we mutate only the elites first?
        idx1, idx2 = sample(1:settings.usersettings.populationsize, 2)
        network1 = population[idx1]
        network2 = population[idx2]
        if network1.fitness < network2.fitness
            push!(newpopulation, deepcopy(network1))
        else
            push!(newpopulation, deepcopy(network2))
        end
    end
    return newpopulation
end

function select_new_population(settings::Settings, population)
    population = sortbyfitness(population)
    newpopulation = eliteselect(settings, population)
    newpopulation  = tournamentselect(settings, population, newpopulation)
    return newpopulation
end

# This is for debugging mostly
function print_top_fitness(n::Int, population)
    for i in 1:n
        println(population[i].fitness)
    end
    return nothing
end


function evolve(settings::Settings, ng::NetworkGenerator, objfunct::ObjectiveFunction)
    population = generate_network_population(settings, ng)
    for i in 1:settings.usersettings.ngenerations
        population = sortbyfitness(population)
        population = mutate_nonelite_population(settings, ng, population)
        population = evaluate_population_fitness(objfunct, population)
        print_top_fitness(1, population)
        population = select_new_population(settings, population)
    end
    population = sortbyfitness(population)
    return population
end