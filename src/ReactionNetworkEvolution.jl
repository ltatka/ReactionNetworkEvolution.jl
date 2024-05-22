module ReactionNetworkEvolution
export run_evolution

import Dates: now, format

include("evo_utils.jl")
include("process_output.jl")

function evolve_networks(batchnum::Int64, parentdir::String, settings::Settings, objfunct::ObjectiveFunction)

    starttime = format(now(), "YYYYmmdd_HHMMSS")
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

    population = generate_network_population(settings, ng)
    species_by_IDs = initialize_species_by_IDs(population)

    if settings.enable_speciation
        species_by_IDs, DELTA = speciate(species_by_IDs, population, DELTA, settings)
    end

    if settings.track_metadata
        num_unique = length(Set(population))
        push!(tracker["num_unique_networks"], num_unique)
        push!(tracker["num_individuals"], length(population))
        push!(tracker["total_num_species"], length(keys(species_by_IDs)))
        avg, mn, mx = get_diversity_stats(species_by_IDs, settings)
        push!(tracker["avg_species_distances"], avg)
        push!(tracker["min_species_distances"], mn)
        push!(tracker["max_species_distances"], mx)
        push!(tracker["delta"], DELTA)
        push!(tracker["avg_intraspecies_distance"], get_intraspecies_distances(species_by_IDs, settings))
    end
    

    for i in 1:settings.ngenerations
    
        # TODO: save time by also returning the top fitness and network?
        species_by_IDs, total_fitness = evaluate_population_fitness(objfunct, species_by_IDs, settings)
        species_by_IDs, total_offspring = calculate_num_offspring(species_by_IDs, total_fitness, settings, writeoutdir=joinpath(starttime, "stalled_models"))


        bestnetwork, maxfitness = gettopmodel(species_by_IDs)


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
            avg, mn, mx = get_diversity_stats(species_by_IDs, settings)
            push!(tracker["avg_species_distances"], avg)
            push!(tracker["min_species_distances"], mn)
            push!(tracker["max_species_distances"], mx)
            push!(tracker["avg_intraspecies_distance"], get_intraspecies_distances(species_by_IDs, settings))
        end

        if settings.enable_speciation
            species_by_IDs, DELTA = speciate(species_by_IDs, population, DELTA, settings)
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
        push!(tracker["avg_intraspecies_distance"], get_intraspecies_distances(species_by_IDs, settings))
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
        stringtracker = json(tracker)
        open(joinpath(starttime, "datatracker.json"), "w") do f
            write(f, stringtracker)
        end
    end

    writeout_settings(settings, joinpath(starttime, "settings.json"))
end



function run_evolution(;
    ngenerations::Int64=-1,
    ntrials::Int64=-1,
    population_size::Int64=-1,
    pathtosettings::String="",
    outputpath::String="",
    seed::Int64=-1,
    note::String="")

    starttime = format(now(), "yyyy-mm-dd_HHMMSS")
    println("Starting $starttime")

    if ntrials == -1
        ntrials = 100
    end
    if outputpath == ""
        outputpath = joinpath(pwd(), "evolution_output")
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

    # If no path to settings is supplied, use default settings
    if pathtosettings == ""
        pathtosettings = "DEFAULT"
    end

    settings, objectivefunction = read_usersettings(pathtosettings, ngenerations=ngenerations, population_size=population_size, seed=seed, note=note)

    for i in 1:ntrials
        evolve_networks(i, path, settings, objectivefunction)
    end

    if settings.process_output_oscillators
        process_oscillators(outputpath)
    else
        println("Output written to $path")
    end

    println("Done")
end

end #module