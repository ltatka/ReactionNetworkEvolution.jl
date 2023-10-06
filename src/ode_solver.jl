using DifferentialEquations
using Sundials
include("reaction_network.jl")
# include("settings.jl")

struct ObjectiveFunction
    objectivespecies::Vector{String} # Can be more than 1
    objectivedata:: DataFrame 
    time::Vector{Float64}
    indexbyspecies::Dict
end

function get_objectivefunction(settings::UserSettings)
    objectivedata = DataFrame(CSV.File(settings.objectivedatapath))
    time = objectivedata[!, "time"]
    indexbyspecies = Dict()
    for s in settings.objectivespecies
        idx = findfirst(item -> item == s, settings.specieslist)
        indexbyspecies[s] = idx
    end
    return ObjectiveFunction(settings.objectivespecies, objectivedata, time, indexbyspecies)
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

function solve_ode(objfunct, network)
    # 
    #Get time info
    t0 = first(objfunct.time)
    t_end = last(objfunct.time)
    stepsize = objfunct.time[2] - t0 # for now we assume even spacing of time points
    tspan = (t0, t_end)

    u0 = network.initialcondition
    ode_prob = ODEProblem(ode_funct!, u0, tspan, network)
    sol = solve(ode_prob, CVODE_BDF(), saveat=stepsize)
    return sol
end

# function evaluate_fitness(objfunct:: ObjectiveFunction, network::ReactionNetwork)
#     sol = solve_ode(objfunct, network)
#     idx = findfirst(item -> item == objfunct.objectivespecies[1], network.specieslist)
#     fitness = 0.0
#     for (i, row) in enumerate(sol.u)
#         fitness += abs(objfunct.objectivedata[1][i] - row[idx])
#     end
#     return fitness # Or should this also assign the fitness to the network?
# end

function evaluate_fitness(objfunct:: ObjectiveFunction, network::ReactionNetwork)
    sol = solve_ode(objfunct, network)
    fitness = 0.0
    for (i, row) in enumerate(sol.u)
        for s in objfunct.objectivespecies
            idx = objfunct.indexbyspecies[s]
            fitness += abs(objfunct.objectivedata[!, s][i] - row[idx])
        end
    end
    return fitness # Or should this also assign the fitness to the network?
end

function evaluate_population_fitness(objfunct::ObjectiveFunction, population)
    for network in population
        fitness = evaluate_fitness(objfunct, network)
        network.fitness = fitness
    end
    return population
end