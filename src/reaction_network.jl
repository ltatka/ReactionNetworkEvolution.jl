using Random
using Distributions
include("settings.jl")

Random.seed!(2)


mutable struct Reaction
    substrate::Vector{String}
    product::Vector{String}
    rateconstant::Float64
end


# struct ReactionProbabilities
#     uniuni::Float64
#     unibi::Float64
#     biuni::Float64
#     bibi::Float64
# end

mutable struct ReactionNetwork
    specieslist::Vector{String}
    initialcondition::Vector{Float64}
    reactionlist::Vector{Reaction}
    floatingspecies::Vector{String}
    boundaryspecies::Vector{String}
    float_initialcondition::Vector{Float64}
    boundary_initialcondition::Vector{Float64}
    fitness::Float64

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
    
    function NetworkGenerator(specieslist::Vector{String}, initialcondition::Vector{Float64}, numreactions::Int, reactionprobabilities::ReactionProbabilities, rateconstantrange::Vector{Float64})
        new(specieslist, initialcondition, numreactions, reactionprobabilities, rateconstantrange)
    end
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
    reactionlist = Vector{Reaction}()
    for i in 1:ng.numreactions
        push!(reactionlist, generate_random_reaction(ng))
    end
    return reactionlist
end

#TODO: Make this more realistic? Maybe a poisson distribution? Or maybe in real life the user will know what species are boundaries?
#TODO: maybe the number of possible boundary species depends on how many species are in the network?
function choose_boundary_species(ng::NetworkGenerator)
    # Decide how many boundary species there will be, 0 to 3, but rarely 3
    weights = [0.6, 0.2, 0.2]
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
    return ReactionNetwork(ng.specieslist, ng.initialcondition, reactionlist, floatingspecies, boundaryspecies, float_initialcondition, boundary_initialcondition, 10E8)
end

# TODO: What about boundary species?
function convert_antimony(network::ReactionNetwork)
    reactions = ""
    rateconstants = ""
    initialconditions = ""
    for (i, reaction) in enumerate(network.reactionlist)
        rxnstring = ""
        ratelawstring = ""
        if length(reaction.substrate) == 1
            rxnstring *= "$(reaction.substrate[1]) -> "
            ratelawstring *= "k$i*$(reaction.substrate[1])"
        elseif length(reaction.substrate) == 2
            rxnstring *= "$(reaction.substrate[1]) + $(reaction.substrate[2]) -> "
            ratelawstring *= "k$i*$(reaction.substrate[1])*$(reaction.substrate[2])"
        else
            rxnstring *= "-> "
        end
        if length(reaction.product) == 1
            rxnstring *= "$(reaction.product[1]); "
        elseif length(reaction.product) == 2
            rxnstring *= "$(reaction.product[1]) + $(reaction.product[2]); "
        else
            rxnstring *= "; "
        end
        rxnstring *= ratelawstring
        reactions *= "$rxnstring\n"
        rateconstants *= "k$i = $(reaction.rateconstant)\n"
    end
    for (i, s) in enumerate(network.specieslist)
        initialconditions *= "$s = $(network.initialcondition[i])\n"
    end
    if length(network.boundaryspecies) > 0
        boundarystr = "const "
        for b in network.boundaryspecies
            boundarystr *= "$b, "
        end
        #remove the extra ", " from the String
        boundarystr = chop(boundarystr, tail=2)
        initialconditions *= boundarystr
    end
    astr = reactions * rateconstants * initialconditions
    return astr
end


