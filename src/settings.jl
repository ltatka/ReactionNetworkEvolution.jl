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
SEED = 6543210

Random.seed!(SEED)

# global global_innovation_number = 0 #TODO: don't use global vars
# global current_innovation_num_by_reaction = Dict()
# delta = 500 #????



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

# TODO: is this structure necessary? Can we just use a dictionary?
struct Settings
    portionelite::Float64
    reactionprobabilities::ReactionProbabilities
    p_rateconstantmutation::Float64
    rateconstantrange::Vector{Float64}
    percent_rateconstant_change::Float64 # Uniform sampling across this range to change rate constant
    p_picknewrateconstant::Float64 # Probability of picking new rate constant during mutation vs slight change
    populationsize::Int
    ngenerations::Int
    nreactions::Int
    max_offspring_portion::Float64 # The maximum portion of offspring a single species can have
    writeout_threshold::Float64 # Networks with this fitness or better will be saved
    p_crossover::Float64 # Probability of crossover vs mutation
    drop_portion::Float64 # Portion of worst networks to drop in each species
    seed::Float64
    starting_delta::Float64
    delta_step::Float64
    target_num_species::Int64
    use_seed_network::Bool # Start with a seed network
    seed_network_path::String 
    # These must be defined in the JSON
    specieslist::Vector{String} # TODO: Maybe they don't need to define this if it's in the data?
    initialconditions::Vector{Float64}
    objectivedatapath::String
    objectivespecies::Vector{String}
   
end


settings = Dict(
    "portionelite" => .1,
    "reactionprobabilities" => [.2, .3, .3, .2],
    "p_rateconstantmutation" => .6, # Probability of changning rate constant vs reaection
    "rateconstantrange" => [0.1, 50.0],
    "percent_rateconstant_change" => 0.2,
    "p_picknewrateconstant" => 0.15,
    "populationsize" => 100,
    "ngenerations" => 400,
    "nreactions" => 5,
    "max_offspring_portion" => 0.1,
    "writeout_threshold" => 0.0088,
    "p_crossover" => 0.2,
    "drop_portion" => 0.1,
    "seed" => -1,
    "starting_delta" => 0.65,
    "delta_step" => 0.1,
    "target_num_species" => 10,
    "use_seed_network" => false,
    "seed_network_path" => ""
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
    settings["seed"] = SEED
    # Create settings object
    # p_rateconstantmutation = p_rateconstantmutation(settings["p_rateconstantmutation"])
    reactionprobabilities = ReactionProbabilities(settings["reactionprobabilities"])
    usersettings = Settings(settings["portionelite"], 
                   reactionprobabilities,
                   settings["p_rateconstantmutation"],
                   settings["rateconstantrange"],
                   settings["percent_rateconstant_change"],
                   settings["p_picknewrateconstant"],
                   settings["populationsize"],
                   settings["ngenerations"],
                   settings["nreactions"],
                   settings["max_offspring_portion"],
                   settings["writeout_threshold"],
                   settings["p_crossover"],
                   settings["drop_portion"],
                   settings["seed"],
                   settings["starting_delta"],
                   settings["delta_step"],
                   settings["target_num_species"],
                   settings["use_seed_network"],
                   settings["seed_network_path"],
                   specieslist,
                   initialconditions,
                   objectivedatapath,
                   objectivespecies)
    return usersettings   
end

function writeout_settings(settings::Settings, filename::String)
    settingsdict = Dict(key=>getfield(settings, key) for key in fieldnames(Settings))
    stringsettings = JSON.json(settingsdict)
    open(filename, "w") do f
        write(f, stringsettings)
    end
end
