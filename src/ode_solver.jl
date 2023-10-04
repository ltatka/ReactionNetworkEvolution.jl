using DifferentialEquations
using Sundials
include("reaction_network.jl")
include("settings.jl")

struct ObjectiveFunction
    objectivespecies::Vector{String} # Can be more than 1
    objectivedata:: Array{Any} 
    time::Vector{Float64}
end

function ode_funct!(du, u, network::ReactionNetwork, t)
    specieslist = network.specieslist
    
    # Reset du
    for i = 1:length(specieslist)
        du[i] = 0.0
    end

    for reaction in network.reactionlist
        # Get the relevant concentrations
        dspecies = 1 # Is there a case where this would be wrong?
        for s in reaction.substrate
            idx = findfirst(item -> item == s, specieslist)
            dspecies *= u[idx]
        end
        # Multiply by the rate constant to get the rate for *this* reaction
        dspecies *= reaction.rateconstant
        # Subtract this rate for substrates
        for s in reaction.substrate
            idx = findfirst(item -> item == s, specieslist)
            du[idx] -= dspecies
        end
        # Add this rate for products
        for p in reaction.product
            idx = findfirst(item -> item == p, specieslist)
            du[idx] += dspecies
        end
    end
    # for boundary species, reset the rate of change to 0
    for s in network.boundaryspecies
        idx = findfirst(item -> item == s, specieslist)
        du[idx] = 0.0
    end
end

function evaluate_fitness(settings::UserSettings, network::ReactionNetwork)

end