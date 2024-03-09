module networkEv

import CSV
import DataFrames
import JSON
import Dates


# include("ode_solver.jl")
include("evo_utils.jl")
# include("network_cleanup.jl")
# include("settings.jl")


println("Starting $(now())")


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
    
    cwd = pwd()
    # pathtosettings = "/home/hellsbells/Desktop/networkEv/test_files/seed_oscillator.json"
    pathtosettings = joinpath(cwd, "test_files/updownObjFunc.json")

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
        # writeoutnetwork(bestnetwork, "$i", directory=joinpath(starttime, "intermediate_models"))

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




num_batches = 100
starttime = now()
starttime = "$starttime"
# Create a directory outside of the networkEv dir to store output
parent = dirname(pwd())
outputdir = joinpath(parent, "Data")
if !isdir(outputdir)
    mkdir(outputdir)
end

# Create a directory of output of this batch in the data dir
path = joinpath(outputdir, "batch_$starttime")
mkdir(path)
println("Writing output to $path")
for i in 1:num_batches
    main(i, path)
end

print("done")

end