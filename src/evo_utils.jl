# include("reaction_network.jl")
include("settings.jl")
using Distributions
using CSV
using DataFrames


function get_objectivefunction(settings::UserSettings)
    df = DataFrame(CSV.File(settings.objectivedatapath))
    time = df[!, "time"]
    objectivedata = Array{Any}[]
    
    for species in settings.objectivespecies
        println(typeof(df[!, species]))
        push!(objectivedata, df[!, species])
    end
    return ObjectiveFunction(settings.objectivespecies, objectivedata, time)
end



function addreaction(ng::NetworkGenerator, network::ReactionNetwork)
    reaction = get_random_reaction(ng)
    push!(network.reactionlist, reaction)
    return network
end

function deletereaction(network::ReactionNetwork)
    idx = sample(1:length(network.reactionlist), 1)
    deleteat!(network.reactionlist, idx)
    return network
end

function mutaterateconstant(network::ReactionNetwork)
    idx = sample(1:length(network.reactionlist), 1)
    percentchange = rand(Uniform(.8, 1.2))
    network.reactionlist[idx][1].rateconstant *= percentchange
    return network
end

