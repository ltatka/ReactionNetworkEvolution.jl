# include("reaction_network.jl")
include("settings.jl")
include("reaction_network.jl")
using Distributions
using CSV
using DataFrames
using Dates
using Random

using Statistics #This is for calculating the mean distance for debugging
using Combinatorics # This is for looking at distances again for debugging



mutable struct Species
    networks::Vector{ReactionNetwork}
    ID::String
    numoffspring::Int64
    speciesfitness::Float64
    topfitness::Float64
    topnetwork::ReactionNetwork
    numstagnations::Int64
    
    #TODO: Now that we've changed the speciate function, these shouldn't all be necessary

    function Species(network::ReactionNetwork)
        # Start a new species with a single network
        networks = [network]
        ID = randstring(12)
        numoffspring = 0
        speciesfitness = 0
        topnetwork = network
        topfitness = 0
        numstagnations = 0
        return new(networks, ID, numoffspring, speciesfitness, topfitness, topnetwork, numstagnations)
    end

    function Species(network::ReactionNetwork, ID::String)
        # Start a new species with a single network and known ID
        networks = [network]
        numoffspring = 0
        speciesfitness = 0
        topnetwork = network
        topfitness = 0
        numstagnations = 0
        return new(networks, ID, numoffspring, speciesfitness, topfitness, topnetwork, numstagnations)
    end

    function Species(network::ReactionNetwork, ID::String, numstagnations::Int64)
        networks = [network]
        numoffspring = 0
        speciesfitness = 0
        topnetwork = network
        topfitness = 0
        return new(networks, ID, numoffspring, speciesfitness, topfitness, topnetwork, numstagnations)
    end


    function Species(listofnetworks::Vector{ReactionNetwork})
        # Start a new species with a list of networks
        networks = listofnetworks
        ID = randstring(12)
        numoffspring = 0
        speciesfitness = 0
        topnetwork = listofnetworks[1] # Obviously this probably isn't going to be the top network at first
        topfitness = 0
        numstagnations = 0
        return new(networks, ID, numoffspring, speciesfitness, topfitness, topnetwork, numstagnations)
    end
end

#TODO: Is this necessary? I thought it would be more complicated.
function initialize_species_by_IDs(population::Array{ReactionNetwork})
    """
    This function is for creating the first species list. It takes a population of networks and assigns them 
    all to the same species. It returns a Species object."""
    species = Species(population)
    ID = species.ID
    return Dict(ID => species)
end


function get_diversity_stats(species_by_IDs::Dict{String, Species})
    # I want to look at average distance between species and maybe also min and max.
    # I'm going to compare every species to every other species 
    all_species = collect(keys(species_by_IDs))
    combo_indices = collect(combinations(1:length(all_species),2))
    distances = []
    for combo in combo_indices
        network1 = species_by_IDs[all_species[combo[1]]].topnetwork
        network2 = species_by_IDs[all_species[combo[2]]].topnetwork
        d = calculate_distance(network1, network2)
        push!(distances, d)
    end
    return mean(distances), minimum(distances), maximum(distances)
end

function get_diversity_stats_from_list(networklist::Array{ReactionNetwork})
    combo_indices = collect(combinations(1:length(networklist),2))
    distances = []
    duplicates = []
    for combo in combo_indices
        network1 = networklist[combo[1]]
        network2 = networklist[combo[2]]
        d = calculate_distance(network1, network2)
        push!(distances, d)
        if d == 0
            push!(duplicates, combo[2])
        end
    end
    duplicates = union(duplicates)
    return mean(distances), minimum(distances), maximum(distances), duplicates
end


function speciate(species_by_IDs::Dict{String, Species},
    population::Vector{ReactionNetwork},
    DELTA::Float64, 
    TARGET_NUM_SPECIES::Int64,
    SPECIES_MOD_STEP::Float64)
    """
    For each network, compare it to the species in the previous generation by 
    randomly picking a network in that species. If it is close enough, then it is a
    member of that species. Otherwise keep looking. If we have checked every species and it is not 
    close enough to any of them, then create a new species for it.

    If there are 2+ network that are similar to each other but dissimilar to any existing species, they will
    each be assigned to a new species. 
    """
    num_species = length(keys(species_by_IDs))
    if num_species > TARGET_NUM_SPECIES
        DELTA += SPECIES_MOD_STEP
    elseif num_species < TARGET_NUM_SPECIES
        DELTA -= SPECIES_MOD_STEP
    end

    new_species_by_IDs = Dict{String, Species}()
    for network in population
        species_assigned = false
        for speciesID in keys(species_by_IDs)
            network2 = rand(species_by_IDs[speciesID].networks)
            distance = calculate_distance(network, network2)
            if distance <= DELTA
                network.ID = speciesID
                # If we've already assigned a network to this species, just add it to the list 
                # Otherwise create a new entry in the hashmap and stop comparing it to existing species (break)
                if speciesID in keys(new_species_by_IDs)
                    push!(new_species_by_IDs[speciesID].networks, network)
                else
                    species = deepcopy(species_by_IDs[speciesID]) # Deepcopy the existing species to retain its data
                    species.networks = [network] # Replace the copies species networks with the new network
                    new_species_by_IDs[speciesID]  = species
                end
                # Once we've assigned a species, stop looking through the old species for a match
                species_assigned = true
                break
            end
        end
        """
        At this point there are two possible states:
            1. We have looked through the previous species, found a match, assigned the new network to it, 
            and broken the loop to land here. In this case, we just want to move onto the next network
            2. We have looked through ALL the previous species and have not found a match. The network has not bestnetwork
            assigned to a species from the previous generation and we must create a new species for it.
        """
        if !species_assigned
            newspecies = Species(network)
            new_species_by_IDs[newspecies.ID] = newspecies
        end
    end
    return new_species_by_IDs, DELTA
end



function getrandomkey(d::Dict{Vector{Vector{String}}, Reaction})
    i = rand(1:length(keys(d)))
    return collect(keys(d))[i]
end


function addreaction(ng::NetworkGenerator, network::ReactionNetwork)
    # if length(network.reactionlist) >= MAX_NREACTIONS
    #     return network
    # end
    reaction = generate_random_reaction(ng)
    # if reaction.key in keys(current_innovation_num_by_reaction)
    #     reaction.innovationnumber = current_innovation_num_by_reaction[reaction.key]
    # else
    #     reaction.innovationnumber = global_innovation_number
    #     current_innovation_num_by_reaction[reaction.key] = global_innovation_number
    #     global_innovation_number += 1
    # end
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


function mutaterateconstant(network::ReactionNetwork, settings::Settings)
    # 90% chance of perturbing existing value, 10% chance of new random value
    key = getrandomkey(network.reactionlist) # Decide which reaction to change
    p = rand()
    if p < settings.p_picknewrateconstant # Randomly pick new rate constant
        newrateconstant = rand(Uniform(settings.rateconstantrange[1], settings.rateconstantrange[2]))
        network.reactionlist[key].rateconstant = newrateconstant
        percentchange = rand(Uniform(.8, 1.2))
        network.reactionlist[key].rateconstant *= percentchange
    else
        percentchange = rand(Uniform(1-settings.percent_rateconstant_change, 1+settings.percent_rateconstant_change))
        network.reactionlist[key].rateconstant *= percentchange
    end
    return network
end

function mutatenetwork!(settings::Settings, ng::NetworkGenerator, network::ReactionNetwork)
    p = rand()
    if p > settings.p_rateconstantmutation
        network = adddeletereaction(ng, network)
    else
        mutaterateconstant(network, settings)
    end
    return network
end

# function mutate_nonelite_population(settings::Settings, ng::NetworkGenerator, population)
#     nelite = Int(floor(settings.portionelite*length(population)))
#     for i in (nelite+1):length(population)
#         population[i] = mutatenetwork!(settings, ng, population[i])
#     end
#     return population
# end

function generate_network_population(settings::Settings, ng::NetworkGenerator)
    population = Vector{ReactionNetwork}()
    # starternetwork = generate_random_network(ng)
    # for i in 1:settings.populationsize
    #     push!(population, deepcopy(starternetwork))
    # end
    for i in 1:settings.populationsize
        network = generate_random_network(ng)
        push!(population, network)
    end
    # # Replace first network with seed network if applicable
    # if ng.seedmodel != Nothing #!isnothing(ng.seedmodel)
    #     population[1] = ng.seedmodel
    # end
    return population
end

function sortbyfitness!(population::Array{ReactionNetwork}, rev::Bool=true)
    # If rev = true, sorts by descending
    # Use rev = true if large fitness = better
    sort!(population, by=get_fitness, rev=rev)
    return population
end


# function eliteselect(settings::Settings, population)
#     nelite = Int(floor(settings.portionelite*length(population)))
#     newpopulation = []
#     for i in 1:nelite
#         push!(newpopulation, deepcopy(population[i]))
#     end
#     return newpopulation
# end


# function tournamentselect(settings::Settings, population, newpopulation)
#     #For now, this is only going to select networks and put them in the new population. Will mutate them later
#     nelite = Int(floor(settings.portionelite*length(population)))
#     for i in 1:(settings.populationsize - nelite)
#         # This will allow elites to be selected and they might dominate
#         # What if we mutate only the elites first?
#         if i < settings.populationsize/2 #for first half, allow elite selection
#             idx1, idx2 = sample(1:settings.populationsize, 2)
#         else # for second half, select only from non-elites
#             if nelite == 0
#                 nelite = 1
#             end
#             idx1, idx2 = sample(nelite:settings.populationsize, 2)
#         end
#         network1 = population[idx1]
#         network2 = population[idx2]
#         if network1.fitness < network2.fitness
#             push!(newpopulation, deepcopy(network1))
#         else
#             push!(newpopulation, deepcopy(network2))
#         end
#     end
#     return newpopulation
# end

# function select_new_population(settings::Settings, population)
#     population = sortbyfitness(population)
#     newpopulation = eliteselect(settings, population)
#     newpopulation  = tournamentselect(settings, population, newpopulation)
#     return newpopulation
# end

# This is for debugging mostly
function print_top_fitness(n::Int, population::Array{ReactionNetwork})
    for i in 1:n
        println(population[i].fitness)
    end
    return nothing
end

function evaluate_population_fitness(objfunct::ObjectiveFunction, species_by_IDs::Dict{String, Species})
    # Evaluates the entire populations fitness as well as the fintess for each species, per the NEAT algo
    # Assigns each individual its fitness score
    # total_fitness = sum of all fitness scores across the entire population, float
    # fitness_by_species = tracks fitness sum of all individuals within a single species, 
    #   hashmap(key, val)= speciesID, float
    # NOTE: The reason this also returns networks_by_species this hashmap is modified, the fitness score is 
    # added to each network
    total_fitness = 0

    for speciesID in keys(species_by_IDs)
        species = species_by_IDs[speciesID]
        species_fitness = 0
        N = length(species.networks)
        topfitness = 0
        topnetwork = species.networks[1]
        for network in species.networks
            fitness = (evaluate_fitness(objfunct, network))
            # total_fitness += fitness
            species_fitness += fitness
            network.fitness = fitness
            # Check if it is the best fitness
            if fitness >= topfitness
                topfitness = fitness
                topnetwork = network
            end
        end
        total_fitness += species_fitness/N
        oldfitness = species.topfitness
        species.topfitness = topfitness
        species.topnetwork = topnetwork
        species.speciesfitness = species_fitness/N #species fitness is now going to be the average fitness of member species
        # If the previous top fitness is the same as this one, increment the stagnation counter by one. Otherwise, reset it to 0
        if oldfitness == topfitness
            species.numstagnations += 1
        else
            species.numstagnations = 0
        end
    end
    return species_by_IDs, total_fitness
end


 

function calculate_num_offspring(species_by_IDs::Dict{String, Species}, total_fitness::Float64, settings::Settings)
    # Calculates how many offspring each species gets based on its share of the total fitness (sum of all individuals)
    # Returns a hashmap key, val = speciesID, float, number of offspring for that species
    total = settings.populationsize
    total_offspring = 0 #This is how many offspring we calculate here, for debugging

    # NEW: cap number of offspring to 10% of population
    MAX_OFFSPRING = round(settings.max_offspring_portion*total)

    for speciesID in keys(species_by_IDs)
        species = species_by_IDs[speciesID]
        ##If the champion of the species has stagnated for more than 15 generations, it won't be allowed to reproduce
        if species.numstagnations >= 50
            species.numoffspring = 0
            networks = sortbyfitness!(species.networks) #TODO: This might not be necessary
            writeoutnetwork(networks[1], "$(networks[1].ID).txt")
            if networks[1].fitness > 0.0088
                println("high fitness network writeout: $(networks[1].ID)")
            end
            
            # println("Writing out a network")
        else
            # if species.numstagnations != 0
            #     println("calculate_num_offspring, stagnation count is $(species.numstagnations)")
            # end
            portion_offspring = species.speciesfitness/total_fitness  # The portion of the next generation this species gets to produce
            numoffspring = round(portion_offspring * total) # The number of offspring this species gets to produce

            # Cap number of offspring
            if numoffspring > MAX_OFFSPRING
                numoffspring = MAX_OFFSPRING
            end
            total_offspring += numoffspring


            species.numoffspring = numoffspring # The number of offspring this species gets to produce
        end
    end
    # println("total offspring: $total_offspring")
    return species_by_IDs
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

function reproduce_networks(species_by_IDs, settings::Settings,
    ng::NetworkGenerator,
    objfunct::ObjectiveFunction;
    generation::Int64=0)
    
    newpopulation = Vector{ReactionNetwork}()

    # println("there are $(length(keys(networks_by_species)))species")

    for speciesID in keys(species_by_IDs)


        species = species_by_IDs[speciesID]
        networks = sortbyfitness!(species.networks)
        # println("species $species top fitness is $(networks[1].fitness)")
        totaloffspringadded = 0
        totaloffspring = species.numoffspring
        # If species is only allowed one offspring, take the most fit individual, mutate it, and either pass on the
        # mutated network or the original, whichever is more fit
        if totaloffspring == 1
            network = networks[1]
            newnetwork = deepcopy(network)
            newnetwork = mutatenetwork!(settings, ng, newnetwork)
            newfitness = evaluate_fitness(objfunct, newnetwork)
            if newfitness > network.fitness
                push!(newpopulation, newnetwork)
            else
                push!(newpopulation, network)
            end
        elseif totaloffspring > 1 # Total offspring greater than 1
            # Directly copy the best network if there are five or more individuals in the species
            if length(networks) >= 5
                push!(newpopulation, deepcopy(networks[1]))
                totaloffspringadded += 1
            end
            # Get rid of the worst networks in the species
            #TODO: Take a look at this if stuff doesn't seem to be working:
            num_to_remove = Int64(floor(length(networks)*settings.drop_portion))
            networks = networks[totaloffspringadded+1:end - num_to_remove] # If we copied over an elite network already, skip it
            # For the rest of the new population:
            offspring_to_add = totaloffspring - totaloffspringadded
            for i in 1:offspring_to_add
                idx1, idx2 = rand(1:length(networks), 2) #TODO: this does NOT prevent same network from being selected twice, which might be good if we only have one more slot to fill
                network = networks[idx1]
                # Decide to mutate it or cross it over
                p = rand()
                if p < settings.p_crossover # crossover with another random network (TODO: for now idc if it crosses over with itself or the elites)
                    network2 = networks[idx2]
                    network2 = networks[idx2]
                    newnetwork = crossover(network, network2)
                    push!(newpopulation, newnetwork)
                else 
                    newnetwork = deepcopy(network)
                    mutatenetwork!(settings, ng, newnetwork)
                    push!(newpopulation, newnetwork)
                end
            end
        end
    end
    return newpopulation
end


function calculate_distance(network1::ReactionNetwork, network2::ReactionNetwork)
    # If every single reaction is different, then the distance will be 1
    W = 0
    num_diff = 0
    c1 = 1 # Value from paper
    c2 = 0.4 # Value from paper (as C3 in the paper) 
    if length(network1.reactionlist) >= length(network2.reactionlist)
        largernetwork = network1
        smallernetwork = network2
    else
        largernetwork = network2
        smallernetwork = network1
    end

    N = length(largernetwork.reactionlist)

    for key in keys(largernetwork.reactionlist)
        if key ∉ keys(smallernetwork.reactionlist)
            num_diff += 1
        end
    end
    
    return num_diff/N

    # for key in keys(network1.reactionlist)
    #     if key in keys(network2.reactionlist)
    #         W += abs(network2.reactionlist[key].rateconstant - network1.reactionlist[key].rateconstant)
    #     else
    #         num_diff += 1
    #     end
    # end

    # return c1*(num_diff)/N + c2*W/N
end

function crossover(network1::ReactionNetwork, network2::ReactionNetwork)
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

function writeoutnetwork(network::ReactionNetwork, filename::String; writeout_threshold::Float64=0.0088, directory::String="stalled_models")
    astr = convert_to_antimony(network)
    astr *= "\n#fitness: $(network.fitness)"

    if !isdir(directory)
        mkdir(directory)
    end
    
    if network.fitness > 0.88
        filename = "highfitness_" * filename
    end
    path = joinpath(directory, filename)

    open(path, "w") do file
        write(file, astr)
    close(file)
    end
end



# function evolve(settings::Settings, ng::NetworkGenerator, objfunct::ObjectiveFunction; writeout=true)
#     # Generate a population consisting of single random network
#     population = generate_network_population(settings, ng)
#     # Make a dictionary of the networks by their species ID
#     networks_by_species = Dict(population[1].ID => population)
#     for i in 1:settings.ngenerations
#         population = sortbyfitness(population)
#         # originset = Set()
#         # for model in population
#         #     push!(originset, model.ID)
#         # end
#         # if length(originset) < 4
#         #     return nothing
#         # end
#         # push!(NUM_VARIANTS, length(originset))
#         # if i%10 == 0
#         #     print_top_fitness(1, population)
#         #     println("set size: $(length(originset))")
#         # end
#         # println(last(population).fitness)
#         population = mutate_nonelite_population(settings, ng, population)
#         population = evaluate_population_fitness(objfunct, population)
#         population = select_new_population(settings, population)
        
#         if writeout
#             fname = "generation_$i.txt"
#             open(fname, "a") do file
#                 for model in population
#                     write(file, "$(model.ID)\n")
#                 end
#             close(file)
#             end
#             fname = "generation_$(i)_fitness"
#             open(fname, "a") do file
#                 for model in population
#                     write(file, "$(model.fitness)\n")
#                 end
#             close(file)
#             end
#         end

#     end
#     population = sortbyfitness(population)
#     return population
# end
