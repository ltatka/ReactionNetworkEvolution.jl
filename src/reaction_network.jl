using Random
using Distributions
using PyCall
include("settings.jl")

Random.seed!(3)


mutable struct Reaction
    substrate::Vector{String}
    product::Vector{String}
    rateconstant::Float64
end

mutable struct ReactionNetwork
    specieslist::Vector{String}
    initialcondition::Vector{Float64}
    reactionlist::Vector{Reaction}
    floatingspecies::Vector{String}
    boundaryspecies::Vector{String}
    float_initialcondition::Vector{Float64}
    boundary_initialcondition::Vector{Float64}
    fitness::Float64
    ID::String # This is for development

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

function get_networkgenerator(settings::Settings)::NetworkGenerator
    specieslist = settings.specieslist
    initialconditions = settings.initialconditions
    nreactions = settings.nreactions
    reactionprobabilities = settings.reactionprobabilities
    rateconstantrange = settings.rateconstantrange
    return NetworkGenerator(specieslist, initialconditions, nreactions, reactionprobabilities, rateconstantrange)
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
    return ReactionNetwork(ng.specieslist, ng.initialcondition, reactionlist, floatingspecies, boundaryspecies, float_initialcondition, boundary_initialcondition, 10E8, ID)
end

function reactants_isequal(reactants1, reactants2)
    if length(reactants1) != length(reactants2)
        return false
    elseif length(reactants1) == 1
        return reactants1 == reactants2
    else
        return (reactants1[1] == reactants2[1] && reactants1[2] == reactants2[2]) ||
            (reactants1[1] == reactants2[2] && reactants1[2] == reactants2[1])

    end
end

function reaction_isnull(reaction)
    # returns true is reaction is pointless, eg. S1 -> S1
    return reactants_isequal(reaction.substrate, reaction.product)
end

function reaction_isequal(reaction1::Reaction, reaction2::Reaction)
    return reactants_isequal(reaction1.substrate, reaction2.substrate) &&
        reactants_isequal(reaction1.product, reaction2.product)    
end 

function removenullreactions(network::ReactionNetwork)
    newreactionlist = []
    for reaction in network.reactionlist
        if !reaction_isnull(reaction)
            push!(newreactionlist, reaction)
        end
    end
    network.reactionlist = newreactionlist
    return network
end

struct ReactionSet
    reactionlist:: Vector{Reaction}
    
    function ReactionSet()
        new([])
    end

end

function pushreaction(reactionset::ReactionSet, reaction; addrateconstants=true)
    for r in reactionset.reactionlist
        if reaction_isequal(r, reaction)
            if addrateconstants
                r.rateconstant += reaction.rateconstant
            end
            return reactionset
        end
    end
    push!(reactionset.reactionlist, reaction)
    return reactionset
end

# TODO: What about boundary species?
function convert_to_antimony(network::ReactionNetwork)
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

function parse_reaction(line::String)
    # Returns reaction object and the name of the rate constant, eg. "k1"
    # Strip spaces
    line = filter(x -> !isspace(x), line)
    if occursin("->", line)
        substrates, products_ratelaw = split(line, "->")
    elseif occursin("=>", line)
        substrates, products_ratelaw = split(line, "=>")
    else
        println("invalid reaction string")
        return nothing # Not a valid reaction, hopefully this line is never called
    end
    # substrates
    s = split("+", substrates)
    products, ratelaw = split(";", products_ratelaw) #TODO: edit this if using more rate laws
    p = split("+", products)
    # Get name of rate parameter
    for subs in s
        ratelaw = replace(ratelaw, s=> "")
    end
    ratelaw = replace(ratelaw, '*'=> "")
    return Reaction(s, p, 0.0), ratelaw
end


function separate_antimony_elements(antstr::String)
    antlines = split(antstr, '\n')
    popfirst!(antlines) # get rid of "// Created by ..." line
    #TODO: make this less sketch?
    all_arrays = [[],[],[],[], []]
    index = 0
    for line in antlines
        if line == ""
            continue
        elseif startswith(line, "//")
            index += 1
            continue
        else
            push!(all_arrays[index], line)
        end
    end
    return all_arrays

end

function process_species_lines(specieslines)
    specieslines = split(specieslines[1], " ")
    specieslist = []
    boundary = []
    floating = []
    for str in specieslines[2:end] # first item will be "species"
        str = replace(str, "," => "")
        str = replace(str, ";" => "")
        if startswith(str, "\$")
            str = str[2:end]
            push!(boundary, str)
        else
            push!(floating, str)
        end
        push!(specieslist, str)
    end
    return specieslist, boundary, floating
end

function process_reaction_lines(reactionlines)
    #TODO: handle stoichiometric coefficients in antimony model?
    reactions_by_ratesymbol = Dict()
    for line in reactionlines
        line = replace(line, "=>" =>"->") # Make sure they all have same arrow
        line = replace(line," "=>"")
        line = replace(line, "\$" => "")
        line = split(line, ":")[2] # Remove reaction name
        line, ratesymbol, _ = split(line, ";") # Remove rate law (for now)
        # Substrates and products
        substrates, products = split(line, "->")
        substrates = split(substrates, "+")
        products = split(products, "+")
        # Get the symbol for the rate constant
        ratesymbol = replace(ratesymbol,"*" => "")
        for s in substrates
            ratesymbol = replace(ratesymbol, s=>"")
        end
        reactions_by_ratesymbol[ratesymbol] = Reaction(substrates, products, 0.0) # We'll deal with rate constant later...?
    end
    return reactions_by_ratesymbol
end

function process_initialcondition_lines(iclines)
    initial_conditions = []
    for line in iclines
        line = replace(line, ";"=>"")
        line = replace(line, " "=>"")
        _, value = split(line, "=")
        push!(initial_conditions, parse(Float64, value))
    end
    return initial_conditions
end

function process_rateconstant_lines(rateconstantlines)
    rate_by_ratesymbol = Dict()
    for line in rateconstantlines
        line = replace(line, ";"=>"")
        line = replace(line, " "=>"")
        symbol, value = split(line, "=")
        rate_by_ratesymbol[symbol] = parse(Float64, value)
    end
    return rate_by_ratesymbol
end

function assign_rates_to_reaction(reactions_by_ratesymbol, rates_by_ratesymbol)
    reactionlist = []
    for k in keys(reactions_by_ratesymbol)
        reactions_by_ratesymbol[k].rateconstant = rates_by_ratesymbol[k]
        push!(reactionlist, reactions_by_ratesymbol[k])
    end
    return reactionlist
end

function get_initialcondition_values(specieslist, initialconditions, sublist)
    sublist_initialconditions = []
    for s in sublist
        idx = findfirst(item -> item == s, specieslist)
        push!(sublist_initialconditions, initialconditions[idx])
    end
    return sublist_initialconditions
end


function convert_from_antimony(filepath::String)
    # Open antimony file and parse contents
    te = pyimport("tellurium")
    rawantstr = read(filepath, String)
    println(rawantstr) #This is the non-standardized antimony string
    r = te.loada(rawantstr)
    antstring = r.getCurrentAntimony() # This is the standardized antimony string
    # Break up lines by category. Eg. species, reactions, etc
    component_array = separate_antimony_elements(antstring)
    # Process each category
    # Get species lists and initial concenetration    
    specieslist, boundaryspecies, floatingspecies = process_species_lines(component_array[1])
    all_initialconditions = process_initialcondition_lines(component_array[3])
    floating_initialcondtions = get_initialcondition_values(specieslist, all_initialconditions, floatingspecies)
    boundary_initialconditions = get_initialcondition_values(specieslist, all_initialconditions, boundaryspecies)
    # Get reactions and rate constants
    reactions_by_ratesymbol = process_reaction_lines(component_array[2])
    rates_by_ratesymbol = process_rateconstant_lines(component_array[4])
    reactionlist = assign_rates_to_reaction(reactions_by_ratesymbol, rates_by_ratesymbol)
    # Construct the ReactionNetwork
    network = ReactionNetwork(specieslist, all_initialconditions, reactionlist, floatingspecies,
        boundaryspecies, floating_initialcondtions, boundary_initialconditions, 10E8, "seedmodel")
    return network
end
    




# function convert_from_antimony(filepath::String)
#     # TODO: make this suck less. So many assumptions about how model will be written.
#     antstring = readlines(filepath)
#     reactionlist = []
#     reactions_by_ratesymbol = {}
#     conc_by_floating = {}
#     conc_by_boundary = {}
#     for line in antstring
#         if startswith(line, "#") || startswith(line, "")
#             continue
#         end
#         # Strip spaces
#         line = filter(x -> !isspace(x), line)
#         # If it's a reaction
#         if occursin("->", line) || occursin("=>", line)
#             if occursin("->", line)
#                 substrates, products_ratelaw = split(line, "->")
#             elseif occursin("=>", line)
#                 substrates, products_ratelaw = split(line, "=>")
#             else
#                 println("invalid reaction string")
#                 return nothing # Not a valid reaction, hopefully this line is never called
#             end
#             # substrates
#             s = split("+", substrates)
#             for sub in s
#                 if startswith(sub, "\$")
#                     conc_by_boundary[chop(sub, head =1, tail =0)] = 0.0
#                 else
#                    conc_by_floating[sub] = 0.0
#                 end
#             end
#             # Products      
#             products, ratelaw = split(";", products_ratelaw) #TODO: edit this if using more rate laws
#             p = split("+", products)
#             for prod in p
#                 if startswith(prod, "\$")
#                     push!(boundaryspeciesn, chop(prod, head =1, tail =0))
#                 else
#                     push!(conc_by_floating, prod)
#                 end
#             end
#             # Get name of rate parameter
#             for subs in s
#                 ratelaw = replace(ratelaw, s=> "")
#             end
#             ratelaw = replace(ratelaw, '*'=> "")
#             r = Reaction(s, p, 0.0)
#             reactions_by_rate_symbol[ratelaw] = r 
#         elseif startswith(line, "k")
#             # For now, I'm going to assume that rate constants start with k and are listed after reactions
#             k, value = split(line, "=")
#             value = parse(Float64, value)
#             reactions_by_ratesymbol[k].rateconstant = value
#         elseif occursin("=", line)
#             # Assuming that if it's not a reaction, not a rate constant assignemnt, but has an equal sign, 
#             # then it must be an initial condition assignemnt
#             species, conc = split(lin, "=")
#             if occursin("\$", species)
#                 species = chop(species, head =1, tail =0)
#             end
#             if species in keys(conc_by_boundary)
#                 conc_by_boundary[species] = parse(Float64, conc)
#             elseif species in keys(conc_by_floating)
#                 conc_by_floating[species] = parse(Float64, conc)
#             else
#                 println("Species $species not found in floats or boundaries")
#                 exit()
#             end
#         else
#             println("Line '$line' not recognized")
#         end
#     end
#     # Putting it all together
#     for k in keys(reactions_by_ratesymbol)
#         push!(reactionlist, reactions_by_ratesymbol[k])
#     end
#     initialconditions = []
#     boundaryspecies = []
#     floatingspecies = []
#     float_initialcondition = []
#     boundary_initialcondition = []
#     specieslist = []
#     for s in keys(conc_by_boundary)
#         conc = conc_by_boundary[s]
#         push!(boundaryspecies, s)
#         push!(specieslist, s)
#         push!(boundary_initialcondition, conc)
#         push!(initialconditions, conc)
#     end
#     for s in keys(conc_by_floating)
#         conc = conc_by_floating[s]
#         push!(floatingspecies, s)
#         push!(specieslist, s)
#         push!(float_initialcondition, conc)
#         push!(initialconditions, conc)
#     end
#     network = ReactionNetwork(specieslist,initialconditions, reactionlist, floatingspecies, boundaryspecies,
#         float_initialcondition, boundary_initialcondition, 10E8)
#     return network
# end

