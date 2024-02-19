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
using Random
import JSON

DEFAULT_FITNESS = 0
MAX_TIME = 10000
MIN_NREACTIONS = 3
MAX_NREACTIONS = 100
SEED = 3

global global_innovation_number = 0 #TODO: don't use global vars
global current_innovation_num_by_reaction = Dict()
delta = 500 #????

Random.seed!(SEED)

struct ReactionProbabilities
    uniuni::Float64
    unibi::Float64
    biuni::Float64
    bibi::Float64

    function ReactionProbabilities(p::Vector{Any})
        new(p[1], p[2], p[3], p[4])
    end

    function ReactionProbabilities(p1, p2, p3, p4)
        new(p1, p2, p3, p4)
    end
end

struct MutationProbabilities
    rateconstant::Float64
    adddeletereaction::Float64

    function MutationProbabilities(p::Vector{Float64})
        new(p[1], p[2])
    end

    function MutationProbabilities(p1::Float64, p2::Float64)
        new(p1, p2)
    end
end


struct Settings
    portionelite::Float64
    reactionprobabilities::ReactionProbabilities
    mutationprobabilities::MutationProbabilities
    specieslist::Vector{String} # TODO: Maybe they don't need to define this if it's in the data?
    initialconditions::Vector{Float64}
    objectivedatapath::String
    objectivespecies::Vector{String}
    populationsize::Int
    ngenerations::Int
    nreactions::Int
    rateconstantrange::Vector{Float64}
end


settings = Dict(
    "populationsize" => 100,
    "ngenerations" => 400,
    "reactionprobabilities" => [.2, .3, .3, .2],
    "mutationprobabilities" => [.6, .4],
    "portionelite" => .1,
    "nreactions" => 5,
    "rateconstantrange" => [0.1, 3.0],
)

function read_usersettings(path::String)
    j = JSON.parsefile(path)
    # Parse required settings
    specieslist = j["specieslist"]
    initialconditions = j["initialconditions"]
    objectivedatapath = j["objectivedatapath"]
    objectivespecies = j["objectivespecies"]
    # Check for any optional args, use defaults if none 
    # If there are user specified args, replace value 
    for k in keys(settings)
        if k in keys(j)
            settings[k] = j[k]
        end
    end
    # Create settings object
    mutationprobabilities = MutationProbabilities(settings["mutationprobabilities"])
    reactionprobabilities = ReactionProbabilities(settings["reactionprobabilities"])
    usersettings = Settings(settings["portionelite"], reactionprobabilities, mutationprobabilities,
        specieslist,initialconditions, objectivedatapath, objectivespecies, settings["populationsize"],
        settings["ngenerations"], settings["nreactions"], settings["rateconstantrange"])
    return usersettings   
end