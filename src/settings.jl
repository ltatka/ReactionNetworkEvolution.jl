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
using JSON
using DataFrames

DEFAULT_FITNESS = 0
MIN_NREACTIONS = 3
MAX_NREACTIONS = 100


struct ReactionProbabilities
    uniuni::Float64
    unibi::Float64
    biuni::Float64
    bibi::Float64

    function ReactionProbabilities(p::Vector{Float64})
        new(p[1], p[2], p[3], p[4])
    end

    function ReactionProbabilities(p1::Float64, p2::Float64, p3::Float64, p4::Float64)
        new(p1, p2, p3, p4)
    end
end


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
    p_crossover::Float64 # Probability of p_crossover
    p_mutation::Float64 # Probability of mutation (does not need to sum to 1 with p_crossover)
    drop_portion::Float64 # Portion of worst networks to drop in each species
    seed::Float64
    starting_delta::Float64
    delta_step::Float64
    target_num_species::Int64
    use_seed_network::Bool # Start with a seed network
    seed_network_path::String
    tournamentselect::Bool
    # These must be defined in the JSON
    specieslist::Vector{String} # TODO: Maybe they don't need to define this if it's in the data?
    initialconditions::Vector{Float64}
    objectivedatapath::String
    #objectivespecies::Vector{String}
   
end

struct ObjectiveFunction
    #objectivespecies::Vector{String} # Can be more than 1
    objectivedata:: Any# DataFrame #TODO: Change this back?
    time::Vector{Float64}
    #indexbyspecies::Dict
end

function get_objectivefunction(settings::Settings)
    if settings.objectivedatapath != "DEFAULT"
        objectivedataframe = DataFrame(CSV.File(settings.objectivedatapath))
        objectivedata = objectivedataframe[!, "S"]
        time = objectivedataframe[!, "time"]
    else
        time = collect(range(0, 1.25, length=11))
        objectivedata = [5.0, 30.0, 5.0, 30.0, 5.0, 30.0, 5.0, 30.0, 5.0, 30.0, 5.0]
    end
    return ObjectiveFunction(objectivedata, time)
end


settings = Dict(
    "portionelite" => 0.1,
    "reactionprobabilities" => [.1, .4, .4, .1],
    "p_rateconstantmutation" => .6, # Probability of changning rate constant vs reaection
    "rateconstantrange" => [0.1, 50.0],
    "percent_rateconstant_change" => 0.2,
    "p_picknewrateconstant" => 0.15,
    "populationsize" => 100,
    "ngenerations" => 800,#400,
    "nreactions" => 5,
    "max_offspring_portion" => 0.1,
    "writeout_threshold" => 0.0088,
    "p_crossover" => 0.75,
    "p_mutation" => 0.75,
    "drop_portion" => 0.1,
    "seed" => -1,
    "starting_delta" => 0.65,
    "delta_step" => 0.1,
    "target_num_species" => 10,
    "use_seed_network" => false,
    "seed_network_path" => "",
    "tournamentselect" => false,
    "specieslist" => ["S0", "S1", "S2"],
    "initialconditions" => [1.0, 5.0, 9.0],
    "objectivedatapath" => "DEFAULT"
)

function read_usersettings(path::String; ngenerations::Int64=-1, populationsize::Int64=-1, seed::Int64=-1)
    # If a path to settings is supplied:
    if path != :"DEFAULT"
        println("Reading settings from $path")
        j = JSON.parsefile(path)
        # # Parse required settings
        # specieslist = j["specieslist"]
        # initialconditions = j["initialconditions"]
        # objectivedatapath = j["objectivedatapath"]

        # Check for any optional args, use defaults if none 
        # If there are user specified args, replace value 
        for k in keys(settings)
            if k in keys(j)
                settings[k] = j[k]
            end
        end
    else
        println("Using default settings")
    end
    # If seed is given as an optional arg, use it and save it to settings. 
    # This will take precedence over any seed specified in the settings file
    # If no random seed is given, check if one is specified in settings, if not, pick randomly and save it
    if seed != -1
        settings["seed"] = seed
    else
        if settings["seed"] == -1
            seed = rand(0:1000000)
            settings["seed"] = seed
        else
            seed = settings["seed"]
        end
    end
    # Set the random seed
    Random.seed!(seed)

    # Create settings object
    reactionprobabilities = ReactionProbabilities(settings["reactionprobabilities"])
    if ngenerations == -1
        ngenerations = settings["ngenerations"]
    end
    if populationsize == -1
        populationsize = settings["populationsize"]
    end
    
    usersettings = Settings(settings["portionelite"], 
                   reactionprobabilities,
                   settings["p_rateconstantmutation"],
                   settings["rateconstantrange"],
                   settings["percent_rateconstant_change"],
                   settings["p_picknewrateconstant"],
                   populationsize,
                   ngenerations,
                   settings["nreactions"],
                   settings["max_offspring_portion"],
                   settings["writeout_threshold"],
                   settings["p_crossover"],
                   settings["p_mutation"],
                   settings["drop_portion"],
                   settings["seed"],
                   settings["starting_delta"],
                   settings["delta_step"],
                   settings["target_num_species"],
                   settings["use_seed_network"],
                   settings["seed_network_path"],
                   settings["tournamentselect"],
                   settings["specieslist"],
                   settings["initialconditions"],
                   settings["objectivedatapath"]
                   )
    return usersettings   
end

function writeout_settings(settings::Settings, filename::String)
    settingsdict = Dict(key=>getfield(settings, key) for key in fieldnames(Settings))
    stringsettings = JSON.json(settingsdict)
    open(filename, "w") do f
        write(f, stringsettings)
    end
end
