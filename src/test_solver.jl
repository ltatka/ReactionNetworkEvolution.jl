using Plots
using CSV
using DataFrames
using JSON
include("ode_solver.jl")
include("evo_utils.jl")
include("network_cleanup.jl")
# include("settings.jl")

using Dates
println("Starting $(now())")

function gettopmodel(species_by_IDs::Dict{String, Species})
    maxfitness = 0
    topnetwork = nothing
    topspecies = nothing
    for speciesID in keys(species_by_IDs)
        species = species_by_IDs[speciesID]
        if species.topfitness > maxfitness
            maxfitness = species.topfitness
            topnetwork = species.topnetwork
            topspecies = species
        end
    end
    return topnetwork, maxfitness
end

function main(batchnum::Int64, parentdir::String)
    starttime = now()
    starttime = "$starttime"
    starttime = joinpath(parentdir, starttime)
    mkdir(starttime)
    tracker = Dict{String, Any}(
        "batch_num" => batchnum,
        "top_individual_fitness" => Vector{Float64}(),
        "total_num_species" => Vector{Int64}(),
        "avg_species_distances" => Vector{Float64}(),
        "min_species_distances" => Vector{Float64}(),
        "max_species_distances" => Vector{Float64}()
    )
    
    # pathtosettings = "/home/hellsbells/Desktop/networkEv/test_files/seed_oscillator.json"
    pathtosettings = "/home/hellsbells/Desktop/networkEv/test_files/updownObjFunc.json"

    settings = read_usersettings(pathtosettings)

    objfunct = get_objectivefunction(settings)
    ng = get_networkgenerator(settings)

    DELTA = settings.starting_delta
    TARGET_NUM_SPECIES = settings.target_num_species
    SPECIES_MOD_STEP = settings.delta_step

    population = generate_network_population(settings, ng)

    species_by_IDs = initialize_species_by_IDs(population)
    species_by_IDs, DELTA = speciate(species_by_IDs, population, DELTA, TARGET_NUM_SPECIES, SPECIES_MOD_STEP)

    species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs)
    species_by_IDs = calculate_num_offspring(species_by_IDs, total_fitness, settings)

    println(settings.ngenerations)

    for i in 1:settings.ngenerations

        species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs)
        species_by_IDs = calculate_num_offspring(species_by_IDs, total_fitness, settings, writeoutdir=joinpath(starttime, "stalled_models"))
        
        L = length(keys(species_by_IDs))
        avg, mn, mx = get_diversity_stats(species_by_IDs)
        push!(tracker["total_num_species"], L)
        push!(tracker["avg_species_distances"], avg)
        push!(tracker["min_species_distances"], mn)
        push!(tracker["max_species_distances"], mx)
        
        bestnetwork, maxfitness = gettopmodel(species_by_IDs)
        writeoutnetwork(bestnetwork, "$i", directory=joinpath(starttime, "intermediate_models"))

        push!(tracker["top_individual_fitness"], maxfitness)
        
        if maxfitness > 0.02
            """It is possible to have oscillators with a lower fitness than this, 
            but seems that any network with this fitness or higher is certainly an oscillator"""
            break 
        end

        population = reproduce_networks(species_by_IDs, settings, ng, objfunct, generation = i)
        species_by_IDs, DELTA = speciate(species_by_IDs, population, DELTA, TARGET_NUM_SPECIES, SPECIES_MOD_STEP)

    end
    
    species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs)
    bestnetwork, maxfitness = gettopmodel(species_by_IDs)

    writeoutnetwork(bestnetwork, "bestmodel_$(bestnetwork.ID)", directory=joinpath(starttime,"final_models"))
    
    for ID in keys(species_by_IDs)
        species = species_by_IDs[ID]
        topnetwork = species.topnetwork
        if topnetwork.ID != bestnetwork.ID
            writeoutnetwork(topnetwork, "$(species.ID)", directory=joinpath(starttime, "final_models"))
        end
    end

    stringtracker = JSON.json(tracker)
    open(joinpath(starttime, "datatracker.json"), "w") do f
        write(f, stringtracker)
    end
    writeout_settings(settings, joinpath(starttime, "settings.json"))
end



num_batches = 1

for i in 1:num_batches
    starttime = now()
    starttime = "$starttime"
    path = joinpath("/home/hellsbells/Desktop/Data/", "batch_$starttime")
    mkdir(path)
    main(i, path)
end



print("done")