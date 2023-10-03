using DifferentialEquations
using Sundials


function get_index(species::String, specieslist::Vector{String})
    # Get the index of a species in a list of species, returns -1 if not in list
    for i in 1:length(specieslist)
        if species == specieslist[i]
            return i
        end
    end
    return  -1
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
            idx = get_index(s, specieslist)
            dspecies *= u[idx]
        end
        # Multiply by the rate constant to get the rate for *this* reaction
        dspecies *= reaction.rateconstant
        # Subtract this rate for substrates
        for s in reaction.substrate
            idx = get_index(s, specieslist)
            du[idx] -= dspecies
        end
        # Add this rate for products
        for p in reaction.product
            idx = get_index(p, specieslist)
            du[idx] += dspecies
        end
    end
    # for boundary species, reset the rate of change to 0
    for s in network.boundaryspecies
        idx = get_index(s, specieslist)
        du[idx] = 0.0
    end
end