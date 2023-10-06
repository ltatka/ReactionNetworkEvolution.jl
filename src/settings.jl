# Things that can have default settings:
# Reaction probabilities
# Rate constant ranges
# Max num generations
# Fitness threshold
# How much a rate constant can change during mutation

# Things that the user has to specify
# Names of species in the network
# INITIAL CONDITIONS
    # Definitely will have this (or close) for the objective species
    # Should initial conditions of other species be things that can change?
# Time series data (doesn't have to be for all species, can(should?) just be for objective function)
# FUTURE: known reactions, possibly their certainty too

# Things we can caclulate/get
# Time ranges, step size, number of data points  
# 

struct UserSettings
    specieslist::Vector{String} # Maybe they don't need to define this if it's in the data?
    initialconditions::Vector{Float64}
    objectivedatapath::String
    objectivespecies::Vector{String}
    populationsize::Int
    ngenerations::Int
end

struct ReactionProbabilities
    uniuni::Float64
    unibi::Float64
    biuni::Float64
    bibi::Float64
end

struct MutationProbabilities
    rateconstant::Float64
    adddeletereaction::Float64
end


struct EvolutionSettings
    portionelite::Float64
    reactionprobabilities::ReactionProbabilities
    mutationprobabilities::MutationProbabilities
    usersettings::UserSettings
end
