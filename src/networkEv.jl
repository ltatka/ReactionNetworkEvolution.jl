module networkEv

import CSV
import DataFrames
import JSON
import Dates

include("evo_utils.jl")

function evolve_networks(batchnum::Int64, parentdir::String, settings::Settings)

    starttime = now()
    starttime = "$starttime"
    starttime = joinpath(parentdir, starttime)
    mkdir(starttime)
    tracker = Dict{String, Any}(
        "batch_num" => batchnum,
        "top_individual_fitness" => Vector{Float64}(undef, settings.ngenerations),
        # "total_num_species" => Vector{Int64}(),
        # "avg_species_distances" => Vector{Float64}(),
        # "min_species_distances" => Vector{Float64}(),
        # "max_species_distances" => Vector{Float64}()
    )

    objfunct = get_objectivefunction(settings)
    ng = get_networkgenerator(settings)

    DELTA = settings.starting_delta
    TARGET_NUM_SPECIES = settings.target_num_species
    SPECIES_MOD_STEP = settings.delta_step

    population = generate_network_population(settings, ng)

    species_by_IDs = initialize_species_by_IDs(population)

    if settings.enable_speciation
        species_by_IDs, DELTA = speciate(species_by_IDs, population, DELTA, TARGET_NUM_SPECIES, SPECIES_MOD_STEP)
    end

    for i in 1:settings.ngenerations

        species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs)
        species_by_IDs, total_offspring = calculate_num_offspring(species_by_IDs, total_fitness, settings, writeoutdir=joinpath(starttime, "stalled_models"))

        # L = length(keys(species_by_IDs))
        # avg, mn, mx = get_diversity_stats(species_by_IDs)
        # push!(tracker["total_num_species"], L)
        # push!(tracker["avg_species_distances"], avg)
        # push!(tracker["min_species_distances"], mn)
        # push!(tracker["max_species_distances"], mx)
        
        bestnetwork, maxfitness = gettopmodel(species_by_IDs)
        tracker["top_individual_fitness"][i] =  maxfitness
        
        if maxfitness > 0.05
            """It is possible to have oscillators with a lower fitness than this, 
            but seems that any network with this fitness or higher is certainly an oscillator"""
            break 
        end

        population = reproduce_networks(species_by_IDs, settings, ng, objfunct, total_offspring)

        if settings.enable_speciation
            species_by_IDs, DELTA = speciate(species_by_IDs, population, DELTA, TARGET_NUM_SPECIES, SPECIES_MOD_STEP)
        else
            species_by_IDs = initialize_species_by_IDs(population)
        end
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



function run_evolution(;
    ngenerations::Int64=-1,
    nbatches::Int64=-1,
    populationsize::Int64=-1,
    pathtosettings::String="",
    outputpath::String="",
    seed::Int64=-1)

    starttime = now()
    starttime = "$starttime"
    println("Starting $starttime")

    if nbatches == -1
        nbatches = 100
    end
    if outputpath == ""
        # Create a directory outside of the networkEv dir to store output
        parent = dirname(pwd())
        outputpath = joinpath(parent, "evolution_output")
    end
    # Create a directory of output of this batch in the data dir
    if !isdir(outputpath)
        mkdir(outputpath)
    end
    path = joinpath(outputpath, "batch_$starttime")
    mkdir(path)
    println("Writing output to $path")

    # If no path to settings is supplied, use default settings
    if pathtosettings == ""
        pathtosettings = "DEFAULT"
    end

    settings = read_usersettings(pathtosettings, ngenerations=ngenerations, populationsize=populationsize, seed=seed)

    for i in 1:nbatches
        evolve_networks(i, path, settings)
    end

    print("done")
end


run_evolution(seed=22, nbatches=351)

end #module