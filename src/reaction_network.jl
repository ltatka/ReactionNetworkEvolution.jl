using Distributions
using PyCall

include("settings.jl")


mutable struct Reaction
    substrate::Vector{String}
    product::Vector{String}
    rateconstant::Float64
    innovationnumber:: Int64
    isactive::Bool
    key::Vector{Vector{String}}
    
    # function Reaction(substrate, product, rateconstant, -1, true)
    #     substrate = sort(substrate)
    #     product = sort(product)
    #     return new(substrate, product, rateconstant, -1, true)
    # end

    function Reaction(substrate, product, rateconstant)
        substrate = sort(substrate)
        product = sort(product)
        return new(substrate, product, rateconstant, -1, true, [substrate, product])
    end

end

function get_reaction_key(reaction::Reaction)
    return [substrate, product]
end


mutable struct ReactionNetwork
    specieslist::Vector{String}
    initialcondition::Vector{Float64}
    reactionlist::Vector{Reaction}
    floatingspecies::Vector{String}
    boundaryspecies::Vector{String}
    float_initialcondition::Vector{Float64}
    boundary_initialcondition::Vector{Float64}
     # This is for development
    fitness::Float64
    ID::String

    # function ReactionNetwork(specieslist::Vector{String},initialcondition::Vector{Float64},reactionlist::Vector{Reaction},floatingspecies::Vector{String},boundaryspecies::Vector{String},float_initialcondition::Vector{Float64},boundary_initialcondition::Vector{Float64})
    #     return new(specieslist, initialcondition, floatingspecies, boundaryspecies, float_initialcondition,boundary_initialcondition, 1.0E8)
    # end
end

#TODO: floating and boundary species
struct NetworkGenerator
    specieslist::Vector{String}
    initialcondition::Vector{Float64}
    numreactions::Int
    reactionprobabilities::ReactionProbabilities
    rateconstantrange::Vector{Float64}
    seedmodel#::Union{ReactionNetwork, Nothing}

    
    function NetworkGenerator(specieslist::Vector{String}, initialcondition::Vector{Float64}, numreactions::Int, reactionprobabilities::ReactionProbabilities, rateconstantrange::Vector{Float64}; seed=nothing)
        new(specieslist, initialcondition, numreactions, reactionprobabilities, rateconstantrange, seed)
    end

    
end

function get_networkgenerator(settings::Settings; seed=Nothing)#::NetworkGenerator
    specieslist = settings.specieslist
    initialconditions = settings.initialconditions
    nreactions = settings.nreactions
    reactionprobabilities = settings.reactionprobabilities
    rateconstantrange = settings.rateconstantrange

    return NetworkGenerator(specieslist, initialconditions, nreactions, reactionprobabilities, rateconstantrange, seed=seed)
   
end

function get_fitness(network::ReactionNetwork)
    return network.fitness
end

#TODO: prevent pointless reactions and ensure connectivity
function generate_random_reaction(ng::NetworkGenerator)
    p = rand()
    k = rand(Uniform(ng.rateconstantrange[1], ng.rateconstantrange[2]))
    
    selectedspecies = sample(ng.specieslist,4)
    subs = [selectedspecies[1]]
    prods = [selectedspecies[2]]

    p_unibi = ng.reactionprobabilities.uniuni + ng.reactionprobabilities.unibi
    p_biuni = p_unibi + ng.reactionprobabilities.biuni

    #NOTE - there is nothing to prevent this from creating pointless reactions, eg. A -> A
    if ng.reactionprobabilities.uniuni < p && p <= p_unibi # unibi
        push!(prods, selectedspecies[3]) # Add a product
    elseif p_unibi < p && p <= p_biuni # biuni
        push!(subs, selectedspecies[3]) # Add a substrate
    elseif p_biuni < p # bibi
        push!(subs, selectedspecies[3]) # Add a substrate
        push!(prods, selectedspecies[4]) # Add a product
    end
    reaction = Reaction(subs, prods, k)
    return reaction
end

function generate_reactionlist(ng::NetworkGenerator)
    global global_innovation_number
    reactionlist = Vector{Reaction}()
    for i in 1:ng.numreactions
        reaction = generate_random_reaction(ng)
        if reaction.key in keys(current_innovation_num_by_reaction)
            reaction.innovationnumber = current_innovation_num_by_reaction[reaction.key]
        else
            reaction.innovationnumber = global_innovation_number
            global_innovation_number += 1
            current_innovation_num_by_reaction[reaction.key] = reaction.innovationnumber
        end

        push!(reactionlist, reaction)
    end
    return reactionlist
end

#TODO: Make this more realistic? Maybe a poisson distribution? Or maybe in real life the user will know what species are boundaries?
#TODO: maybe the number of possible boundary species depends on how many species are in the network?
function choose_boundary_species(ng::NetworkGenerator)
    # Decide how many boundary species there will be, 0 to 3, but rarely 3
    weights = [1, 0, 0]
    nboundary = wsample(0:2, weights)
    # Decide which will be the boundary species
    boundaryspecies = sample(ng.specieslist, nboundary, replace=false)
    floatingspecies = [s for s in ng.specieslist if !(s in boundaryspecies)]
    return floatingspecies, boundaryspecies
end

#TODO: Might not need to do this for the floats?
function get_initialconditions(ng::NetworkGenerator, floatingspecies::Vector{String}, boundaryspecies::Vector{String})
    float_initialcondition = Vector{Float64}()
    boundary_initialcondition = Vector{Float64}()
    for species in floatingspecies
        idx = findfirst(item -> item == species, ng.specieslist)
        push!(float_initialcondition, ng.initialcondition[idx])
    end
    for species in boundaryspecies
        idx = findfirst(item -> item == species, ng.specieslist)
        push!(boundary_initialcondition, ng.initialcondition[idx])
    end
    return float_initialcondition, boundary_initialcondition
end

function generate_random_network(ng::NetworkGenerator)
    floatingspecies, boundaryspecies = choose_boundary_species(ng)
    float_initialcondition, boundary_initialcondition = get_initialconditions(ng, floatingspecies, boundaryspecies)
    reactionlist = generate_reactionlist(ng)
    ID = randstring(10)
    return ReactionNetwork(ng.specieslist, ng.initialcondition, reactionlist, floatingspecies, boundaryspecies, float_initialcondition, boundary_initialcondition, DEFAULT_FITNESS, ID)
end
