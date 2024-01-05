# include("reaction_network.jl")
include("settings.jl")
using Distributions
using CSV
using DataFrames
using Dates



NUM_VARIANTS = [] 

function getrandomkey(d::Dict)
    i = rand(1:length(keys(d)))
    return collect(keys(d))[i]
end


function addreaction(ng::NetworkGenerator, network::ReactionNetwork)
    global global_innovation_number
    if length(network.reactionlist) >= MAX_NREACTIONS
        return network
    end
    reaction = generate_random_reaction(ng)
    if reaction.key in keys(current_innovation_num_by_reaction)
        reaction.innovationnumber = current_innovation_num_by_reaction[reaction.key]
    else
        reaction.innovationnumber = global_innovation_number
        current_innovation_num_by_reaction[reaction.key] = global_innovation_number
        global_innovation_number += 1
    end
    if reaction.key in keys(network.reactionlist)
        network.reactionlist[reaction.key].rateconstant += reaction.rateconstant 
    else
        network.reactionlist[reaction.key] = reaction
    end
    return network
end

function deletereaction(network::ReactionNetwork)
    if length(network.reactionlist) <= MIN_NREACTIONS # Maintain a minimum of 2 reactions
        return network
    end
    key = getrandomkey(network.reactionlist)
    #TODO: should reactions that are already inactive be excluded from this?
    network.reactionlist[key].isactive = false
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
    key = getrandomkey(network.reactionlist)
    percentchange = rand(Uniform(.8, 1.2))
    network.reactionlist[key].rateconstant *= percentchange
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
    starternetwork = generate_random_network(ng)
    for i in 1:settings.populationsize
        push!(population, deepcopy(starternetwork))
    end
    # for i in 1:settings.populationsize
    #     network = generate_random_network(ng)
    #     push!(population, network)
    # end
    # # Replace first network with seed network if applicable
    # if ng.seedmodel != Nothing #!isnothing(ng.seedmodel)
    #     population[1] = ng.seedmodel
    # end
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
    for i in 1:(settings.populationsize - nelite)
        # This will allow elites to be selected and they might dominate
        # What if we mutate only the elites first?
        if i < settings.populationsize/2 #for first half, allow elite selection
            idx1, idx2 = sample(1:settings.populationsize, 2)
        else # for second half, select only from non-elites
            if nelite == 0
                nelite = 1
            end
            idx1, idx2 = sample(nelite:settings.populationsize, 2)
        end
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

function evaluate_population_fitness(objfunct::ObjectiveFunction, networks_by_species)
    # Evaluates the entire populations fitness as well as the fintess for each species, per the NEAT algo
    # Assigns each individual its fitness score
    # total_fitness = sum of all fitness scores across the entire population, float
    # fitness_by_species = tracks fitness sum of all individuals within a single species, 
    #   hashmap(key, val)= speciesID, float
    total_fitness = 0
    fitness_by_species = Dict()
    for species in keys(networks_by_species)
        species_fitness = 0
        N = length(networks_by_species[species])
        for network in networks_by_species[species]
            fitness = (1/N)*(1/evaluate_fitness(objfunct, network)) # Invert the fitness to make it larger=better
            total_fitness += fitness
            species_fitness += fitness
            network.fitness = fitness
        end
        fitness_by_species[species] = species_fitness
    end
    return networks_by_species, fitness_by_species, total_fitness
end

function calculate_num_offspring(fitness_by_species, total_fitness, settings::Settings)
    # Calculates how many offspring each species gets based on its share of the total fitness (sum of all individuals)
    # Returns a hashmap key, val = speciesID, float, number of offspring for that species
    total = settings.populationsize
    total_calculated = 0
    numoffspring_by_species = Dict()
    for species in keys(fitness_by_species)
        numoffspring = round(fitnress_by_species[species]/total_fitness)
        total_calculated += numoffspring
        numoffspring_by_species[species] = numoffspring
    end
    if total_calculated != total
        println("calculated $total_calculated not eqault to popsize of $total")
    end
    return numoffspring_by_species
end


function evolve(settings::Settings, ng::NetworkGenerator, objfunct::ObjectiveFunction; writeout=true)
    # Generate a population consisting of single random network
    population = generate_network_population(settings, ng)
    # Make a dictionary of the networks by their species ID
    networks_by_species = Dict(population[1].ID => population)
    for i in 1:settings.ngenerations
        population = sortbyfitness(population)
        # originset = Set()
        # for model in population
        #     push!(originset, model.ID)
        # end
        # if length(originset) < 4
        #     return nothing
        # end
        # push!(NUM_VARIANTS, length(originset))
        # if i%10 == 0
        #     print_top_fitness(1, population)
        #     println("set size: $(length(originset))")
        # end
        # println(last(population).fitness)
        population = mutate_nonelite_population(settings, ng, population)
        population = evaluate_population_fitness(objfunct, population)
        population = select_new_population(settings, population)
        
        if writeout
            fname = "generation_$i.txt"
            open(fname, "a") do file
                for model in population
                    write(file, "$(model.ID)\n")
                end
            close(file)
            end
            fname = "generation_$(i)_fitness"
            open(fname, "a") do file
                for model in population
                    write(file, "$(model.fitness)\n")
                end
            close(file)
            end
        end

    end
    population = sortbyfitness(population)
    return population
end

function speciate(networks_by_species, population, delta)
    #global delta
    new_networks_by_species = Dict()
    for network in population
        species_assigned = false
        for species in keys(networks_by_species)
            network2 = rand(networks_by_species[species])
            distance = calculate_distance(network, network2)
            if distance <= delta
                network.ID = species
                if species in keys(new_networks_by_species)
                    push!(new_networks_by_species[species], network)
                else
                    new_networks_by_species[species] = [network]
                end
                species_assigned = true
                break
            end
        end
        # If we've gone through all existing species and none are close enough
        # to the network, then create a new species#TODO WHat if there are two species
        # that are close to one another, but not any existing species? In this case, they will each be species_assigned
        # a new species number. But I can't think of a way to avoid this...
        if !species_assigned
            new_ID = randstring(10)
            new_networks_by_species[new_ID] = [network]
        end
    end
    return new_networks_by_species
end

function cleanupreactions(network::ReactionNetwork)
    network = removenullreactions(network)
    reactionset = ReactionSet()
    for reaction in network.reactionlist
        reactionset = pushreaction(reactionset, reaction)
    end
    network.reactionlist = reactionset.reactionlist
    return network
end


function calculate_distance(network1, network2)
    W = 0
    num_diff = 0
    c1 = 1
    c2 = 1 
    if length(network1.reactionlist) > length(network2.reactionlist)
        N = length(network1.reactionlist)
    else
        N = length(network2.reactionlist)
    end

    for key in keys(network1.reactionlist)
        if key in keys(network2.reactionlist)
            W += abs(network2.reactionlist[key].rateconstant - network1.reactionlist[key].rateconstant)
        else
            num_diff += 1
        end
    end

    return c1*(num_diff)/N + c2*W/N
end

function crossover(network1, network2)
    newreactiondict = Dict()
    if network1.fitness > network2.fitness
        morefitnetwork = network1
        lessfitnetwork = network2
    elseif network1.fitness < network2.fitness
        morefitnetwork = network2
        lessfitnetwork = network1
    else #If equal fitness, randomly assign roles  
        p = rand()
        if p < 0.5
            morefitnetwork = network1
            lessfitnetwork = network2
        else
            morefitnetwork = network2
            lessfitnetwork = network1
        end
    end
    # We are going to look through the genes in the more fit network. If there is a gene in the 
    # less fit network that is not in the more fit one, we don't care. But if there is an unmatched gene in the more
    # fit network, we want to keep it
    for key in keys(morefitnetwork.reactionlist)
        # If the reaction is in both networks, randomly copy it from either network
        if key in keys(lessfitnetwork.reactionlist)
            p = rand()
            if p < 0.5
                newreaction = morefitnetwork.reactionlist[key]
            else
                newreaction = lessfitnetwork.reactionlist[key]
            end
            
        else # If the reaction is NOT in the less fit network, copy it over
            newreaction = morefitnetwork.reactionlist[key]
        end
        # If the selected reaction is inactive, 25% of it being reactivated
        if !newreaction.isactive
            p = rand()
            if p < 0.25
                newreaction.isactive = true
            end
        end
        newreactiondict[key] = newreaction
    end

    newnetwork = deepcopy(morefitnetwork)
    newnetwork.reactionlist = newreactiondict

    #TODO: should we reset the fitness, keep the old one and then replace it later or doesn't matter?
    return newnetwork

end
