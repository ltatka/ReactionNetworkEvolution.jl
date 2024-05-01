module networkEv

import CSV
import DataFrames
import JSON
import Dates

include("evo_utils.jl")

function evolve_networks(batchnum::Int64, parentdir::String, settings::Settings, objfunct::ObjectiveFunction)

    starttime = now()
    starttime = "$starttime" * randstring(3)  # for race conditions
    starttime = joinpath(parentdir, starttime)
    mkdir(starttime)
    if settings.track_metadata
        tracker = Dict{String, Any}(
            "batch_num" => batchnum,
            "top_individual_fitness" => Vector{Float64}(),
            "num_unique_networks" => Vector{Int64}(),
            "num_individuals" => Vector{Int64}(),
            "num_best_network" => Vector{Int64}(),
            "total_num_species" => Vector{Int64}(),
            "avg_species_distances" => Vector{Float64}(),
            "min_species_distances" => Vector{Float64}(),
            "max_species_distances" => Vector{Float64}(),
            "delta" => Vector{Float64}(),
            "avg_intraspecies_distance" => Vector{Float64}()
        )
    end
    
    ng = get_networkgenerator(settings)

    DELTA = settings.starting_delta
    TARGET_NUM_SPECIES = settings.target_num_species
    SPECIES_MOD_STEP = settings.delta_step

    population = generate_network_population(settings, ng)
    species_by_IDs = initialize_species_by_IDs(population)

    if settings.enable_speciation
        species_by_IDs, DELTA = speciate(species_by_IDs, population, DELTA, TARGET_NUM_SPECIES, SPECIES_MOD_STEP)
    end

    if settings.track_metadata
        num_unique = length(Set(population))
        push!(tracker["num_unique_networks"], num_unique)
        push!(tracker["num_individuals"], length(population))
        push!(tracker["total_num_species"], length(keys(species_by_IDs)))
        avg, mn, mx = get_diversity_stats(species_by_IDs)
        push!(tracker["avg_species_distances"], avg)
        push!(tracker["min_species_distances"], mn)
        push!(tracker["max_species_distances"], mx)
        push!(tracker["delta"], DELTA)
        push!(tracker["avg_intraspecies_distance"], get_intraspecies_distances(species_by_IDs))
    end
    

    for i in 1:settings.ngenerations
    
        # TODO: save time by also returning the top fitness and network?
        species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs, settings)
        species_by_IDs, total_offspring = calculate_num_offspring(species_by_IDs, total_fitness, settings, writeoutdir=joinpath(starttime, "stalled_models"))


        bestnetwork, maxfitness = gettopmodel(species_by_IDs)
        if length(tracker["top_individual_fitness"]) > 0 &&  tracker["top_individual_fitness"][end] - maxfitness > 0.0002
            println("generation $i")
        end

        if settings.track_metadata
            push!(tracker["top_individual_fitness"],maxfitness)
            bestnetwork_count = count_best_networks(bestnetwork, population)
            push!(tracker["num_best_network"], bestnetwork_count)
        end


        if maxfitness > settings.writeout_threshold #0.05
            """It is possible to have oscillators with a lower fitness than this, 
            but seems that any network with this fitness or higher is certainly an oscillator"""
            break 
        end

        population = reproduce_networks(species_by_IDs, settings, ng, objfunct, total_offspring)

        if settings.track_metadata
            num_unique = length(Set(population))
            push!(tracker["num_unique_networks"], num_unique)
            push!(tracker["num_individuals"], length(population))
            push!(tracker["total_num_species"], length(keys(species_by_IDs)))
            avg, mn, mx = get_diversity_stats(species_by_IDs)
            push!(tracker["avg_species_distances"], avg)
            push!(tracker["min_species_distances"], mn)
            push!(tracker["max_species_distances"], mx)
            push!(tracker["avg_intraspecies_distance"], get_intraspecies_distances(species_by_IDs))
        end

        if settings.enable_speciation
            species_by_IDs, DELTA = speciate(species_by_IDs, population, DELTA, TARGET_NUM_SPECIES, SPECIES_MOD_STEP)
        else
            species_by_IDs = initialize_species_by_IDs(population)
        end
    end
    
    species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs, settings)
    bestnetwork, maxfitness = gettopmodel(species_by_IDs)

    if settings.track_metadata
        push!(tracker["top_individual_fitness"],maxfitness)
        bestnetwork_count = count_best_networks(bestnetwork, population)
        push!(tracker["num_best_network"], bestnetwork_count)
        push!(tracker["delta"], DELTA)
        push!(tracker["avg_intraspecies_distance"], get_intraspecies_distances(species_by_IDs))
    end


    writeoutnetwork(bestnetwork, "bestmodel_$(bestnetwork.ID)", directory=joinpath(starttime,"final_models"))
    
    for ID in keys(species_by_IDs)
        species = species_by_IDs[ID]
        topnetwork = species.topnetwork
        if topnetwork.ID != bestnetwork.ID
            writeoutnetwork(topnetwork, "$(species.ID)", directory=joinpath(starttime, "final_models"))
        end
    end

    if settings.track_metadata
        stringtracker = JSON.json(tracker)
        open(joinpath(starttime, "datatracker.json"), "w") do f
            write(f, stringtracker)
        end
    end

    writeout_settings(settings, joinpath(starttime, "settings.json"))
end



function run_evolution(;
    ngenerations::Int64=-1,
    nbatches::Int64=-1,
    populationsize::Int64=-1,
    pathtosettings::String="",
    outputpath::String="",
    seed::Int64=-1,
    note::String="")

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
        # This is for dealing with race conditions in parallel computing
        try
            mkdir(outputpath)
        catch LoadError
        end
    end

    path = joinpath(outputpath, "batch_$starttime" * randstring(3)) # randstring is for race conditions

    if !isdir(path)
        mkdir(path)
    end

    println("Writing output to $path")

    # If no path to settings is supplied, use default settings
    if pathtosettings == ""
        pathtosettings = "DEFAULT"
    end

    settings, objectivefunction = read_usersettings(pathtosettings, ngenerations=ngenerations, populationsize=populationsize, seed=seed, note=note)

    for i in 1:nbatches
        evolve_networks(i, path, settings, objectivefunction)
    end

    print("done")
end

end #module